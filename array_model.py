from sympy import *
import numpy as np
import matplotlib.pyplot as plt
import csv
import os
import pandas as pd
import ltspice
import time
import ast
import math

path_here = os.path.dirname(os.path.realpath(__file__))

params = {"ND": 1e16, #provide number if doping type is constant, 1D array if alternating
          "Doping type": "poisson", #options: constant (provide a number for ND), alternating (provide an array), poisson (provide doping mean and length will be taken as constant using mean length)
          "Previous doping distribution":None, #None or provide address to file from this directory (like length one below)
          "NT": 4e11, #default 5e11

          "Contact type": "point", #point or electrode (if point, Input locations, Voltage inputs and Ground locations are IGNORED)
          "Input locations": {"left":([14],["V1"]), "top":None, "right":None, "bottom":None},
          "Voltage inputs" : {"V1":0.5,"V2":0.5,"V3":0, "V4":0, "V5":0},
          "Ground locations":{"left":0, "top":0, "right":[0], "bottom":None}, 
          # example: (left/right from 0 to n-1 (in range n) and top/bottom from 0 to m+1**fix this (in range m+2) - one more node than grain? 
          # NOTE: FORMAT ([voltage input locations], (corresponding voltage names)) for input; [ground locations] for output; can also use "edge" for inputs/outputs
          # "Input locations": {"left":([0, range(1,4)], ("V1", "V2")), "top":None, "right":None, "bottom":None},
          # "Ground locations":{"left":0, "top":0, "right":[range(1,4)], "bottom":None},
          # Inputs only implemented for left, outputs only implemented for right
          "VA LUT": (-1.3, 1.3), #Voltage range used to generate LUT (only one boundary), should slightly exceed the voltage limit if possible (must not be under unless the GB voltage won't exceed that)
          "VA": (-5, 15), #Voltage range used in simulation

          "Grain size distribution": "exponential", #options: constant, exponential
          "Thickness": 1.8e-7, #180 nm
          "Length": 1.5e-3, #for constant: provide length; for exponential: provide mean length, NOTE vary cutoff below
          "Lower grain size limit": 0.5e-4,
          "Upper grain size limit": 100e-4,
          "Rows of grains (n)": 5,
          "Columns of grains (m)": 5,

          "SimDir": 'Report things Appendix netlist',
          "SimFile": 'Si 5x5 ND 1e16 L 15um T 300K NT 4e11 Point Contacts (low res 10)',
          "LUT_filenum": 'lowres10',
          "Previous Length distribution":None, #None if no previous table, or provide its address from this directory, MAKE SURE L dimensions match desired dimensions
            #eg. "Variable length tests\\length_dist_Test_5x5_constant_doping_exp_length_1.5e-3.txt"

          "Electron mobility": 8000, #This is ignored when choosing Silicon, NOTE: model doesn't take into account temperature dependance of mobility
          "Material": "Si", #"Si" for silicon, otherwise provide electron mobility for other materials
          "er":11.7,
          "Temperature":300,
          "Effective mass":1.09, #Calculates Richardson constant and NC (been using 1.09 for Si)
          
          "LUT resolution":10, #PREVIOUS FOR REPORT GRAPHS: 400
          "Voltage sweep step":0.1}

reduced_keys = ["NT", "VA LUT", "Thickness", "SimDir", "LUT_filenum", "er", "Temperature", "Effective mass", "LUT resolution"]
reduced_params = {key: params.get(key) for key in reduced_keys}
title = f"LUT_filenum = {params['LUT_filenum']}, previous doping distibution = {params['Previous doping distribution']}, previous length distribution = {params['Previous Length distribution']}"

def initial_params(T=params['Temperature'], mn_eff=params['Effective mass']):
    NC = 2.5*(mn_eff)**1.5*(T/300)**1.5*10**19
    params['NC'] =  NC
    AA = 120*mn_eff
    params["Richardson constant (AA)"] = AA

initial_params()

def find_widths(T=params['Temperature'], er=params['er'], ND1=1e16, ND2=1e16, NT=params['NT'], VA=0): #finds depletion width
    q = 1.6e-19 #electron charge
    e0 = 8.85e-14 #permittivity of free space
    ee = e0*er # er permittivity of silicon 11.7
    k = 8.6e-5 #boltzmann constant in eV
    kT = k*T

    x1, x2 = symbols('x1 x2')

    eqn1 = Eq(x1*ND1+x2*ND2, NT) #balanced charges
    eqn2 = Eq(q*ND1*x1**2/(2*ee), q*ND2*x2**2/(2*ee) + kT*ln(ND1/ND2) - VA) #potentials meet at boundary (x=0)
    result = solve([eqn1, eqn2], (x1, x2)) #solve simultaneously
    if result[0][0] >= 0 and result[0][1] >= 0:
        x1s = float(result[0][0])
        x2s = float(result[0][1])
    elif len(result)==1:
        print("ERROR: Only one solution. Voltage limit exceeded.")
        x1s = np.nan
        x2s = np.nan
    elif result[1][0] >= 0 and result[1][1] >= 0:
        x1s = float(result[1][0])
        x2s = float(result[1][1])
    else: 
        x1s = np.nan
        x2s = np.nan
        print("ERROR: Voltage limit exceeded.")
    # print("x1 = "+str(x1s)+" and x2 = "+str(x2s)+"\n")

    return x1s, x2s

def find_barriers(x1, x2, T=params['Temperature'], ND1=1e16, ND2=1e16, NC=params['NC'], er=params["er"]):
    q = 1.6e-19
    e0 = 8.85e-14
    k = 8.6e-5
    kT = k*T
    ee = e0*er

    phi_1 = q*ND1*x1**2/(2*ee) + kT*np.log(NC/ND1) #calculating Ec - Efx
    phi_2 = q*ND2*x2**2/(2*ee) + kT*np.log(NC/ND2)
    # print(f'phi_1={phi_1}, phi_2={phi_2}, diff={phi_1-phi_2}')
    return phi_1, phi_2

def calc_current(phi_1, phi_2, A=1, T=300, AA=120):
    k = 8.6e-5
    kT = k*T
    q = 1.6e-19
    current = A*AA*(T**2)*(np.exp(-phi_1/kT) - np.exp(-phi_2/kT))
    return current

def LUTwriter(SimDir=params['SimDir'], filenum=params['LUT_filenum'], ND1=1e16, ND2=1e16, L1=1e-4, L2=1e-4, NT=params['NT'], Vrange=params['VA LUT'], er=params['er'], T=params['Temperature'], res=params['LUT resolution'], NC=params['NC'], AA=params['Richardson constant (AA)'], parameters=reduced_params, t=params['Thickness']):
    L = min(L1, L2)
    A = L*t
    path = path_here+'\\'+SimDir
    if not os.path.exists(path):
        os.mkdir(path)
    LUT_path = path+'\\Lookup Tables'
    if not os.path.exists(LUT_path):
        os.mkdir(LUT_path)
    VA = np.linspace(Vrange[0], Vrange[1], res) #makes array from -1 to 1 with 100 values (res for resolution)
    I = []
    VA_solutions = []
    for k in range(len(VA)):
        # print("VA="+str(VA[k])+": ")
        x1, x2 = find_widths(VA=VA[k], ND1=ND1, ND2=ND2, NT=NT, T=T, er=er)
        if not np.isnan(x1):
            phi_1, phi_2 = find_barriers(x1, x2, ND1=ND1, ND2=ND2, T=T, er=er, NC=NC)
            # print("phi_1="+str(phi_1)+" phi_2="+str(phi_2))
            I.append(calc_current(phi_1, phi_2, A=A, T=T, AA=AA))
            VA_solutions.append(VA[k])
        else:
            print(f'WARNING: No solution at VA={VA[k]}')
    VA_solutions = np.array(VA_solutions)
    I = np.array(I)
    R = VA_solutions/I
    with open(LUT_path+'\\LUT_'+str(filenum)+"_ND_"+str(f'{ND1:.2e}_{ND2:.2e}')+"_L_"+str(f'{L:.2e}')+'.txt', 'w', newline='') as csvfile:
        for key,value in parameters.items():
            if type(value)==str:
                csvfile.write(f'#{key}:\'{value}\'\n')
            else:
                csvfile.write(f'#{key}:{value}\n')
        Lwriter = csv.writer(csvfile, delimiter=' ')
        Lwriter.writerow(['VA', 'R', 'I'])
        Lwriter.writerow
        for k in range(len(VA_solutions)):
            Lwriter.writerow([str(VA_solutions[k]), str(R[k]), str(I[k])])
    print(str(f'LUT made: {LUT_path}\\LUT_{filenum}_ND_{ND1:.2e}_{ND2:.2e}_L_{L:.2e}.txt'))

def density_states(T=params['Temperature'], mn_eff=params['Effective mass']):
    NC = 2.5*(mn_eff)**1.5*(T/300)**1.5*10**19
    params['NC'] =  NC

def bulk_resistance(W, ND=params['ND'], t=params['Thickness'], L=params['Length'], mu_n=params['Electron mobility'], material=params['Material']):
    q =  1.6e-19
    A = t*W
    if material != "Si":
        mu_mono=mu_n
    else:
        mu_mono=1400/math.sqrt(1+(ND/3e16)*350/(ND/3e16+350))
    # print(f'mu_mono={mu_mono}')
    R = L/(q*A*(ND*mu_mono))
    return R

def doping_array(ND=params['ND'], type=params['Doping type'], n=params['Rows of grains (n)'], m=params['Columns of grains (m)'], L=params['Length'], t=params['Thickness'], SimDir=params['SimDir'], filenum=params['SimFile'], load_prev=params['Previous doping distribution']):
    if load_prev!=None:
        ND_array = np.loadtxt(str(f'{path_here}\\{load_prev}'))
        dimensions = ND_array.shape
        if (dimensions[0]!=n) and (dimensions[1]!=2*m-2):
            print("WARNING: Loaded doping array does not match requested dimensions.")
    elif type == "constant" or type == "Constant":
        ND_array = np.full((n,2*m-2),ND)
        # # DEBUG:
        print(str(f'ND_array = {ND_array}'))
    elif type == "alternating" or type =="Alternating":
        ND_array = np.zeros((n,2*m-2))

        #Generate one row of doping concentrations
        ND_row = [ND[0]]
        grain_counter = 1
        for i in range(1,m-1): #do middle grains (could just do 0,m-2)
            if grain_counter < len(ND):
                ND_row.append(ND[grain_counter])
                ND_row.append(ND[grain_counter])
                grain_counter += 1
            else:
                grain_counter = 0
                ND_row.append(ND[grain_counter])
                ND_row.append(ND[grain_counter])
                grain_counter += 1

        if grain_counter < len(ND):
            ND_row.append(ND[grain_counter])
        else:
            ND_row.append(ND[0])

        print(str(f'DEBUG: ND_row = {ND_row}'))

        for i in range(n):
            ND_array[i] = ND_row
            ND_row.append(ND_row[len(ND_row)-1])
            ND_row.append(ND_row[0])
            del(ND_row[0:2])
            # print(str(f'ND_row = {ND_row}'))
        
        print(str(f'DEBUG: Doping array = {ND_array}'))
    elif type == "poisson" or type == "Poisson": 
        volume = t*L**2
        atoms_per_grain = ND*volume
        atoms_poisson=np.random.poisson(atoms_per_grain, size=(n,m))
        ND_poisson = atoms_poisson/volume
        ND_array=np.zeros((n,2*m-2))
        ND_array[:,0]=ND_poisson[:,0]
        ND_array[:,-1]=ND_poisson[:,-1]
        grain_counter=1
        for i in range(1, 2*m-3, 2):
            ND_array[:,i:i+2]=np.column_stack((ND_poisson[:, grain_counter], ND_poisson[:, grain_counter]))
            grain_counter+=1
        np.savetxt(str(f'{path_here}\\{SimDir}\\poisson_doping_dist_{filenum}.txt'), ND_array)
        print(str(f'Doping array = {ND_array}'))
    return ND_array

def length_array(L_mean = params['Length'], n=params['Rows of grains (n)'], m=params['Columns of grains (m)'], dist_type=params['Grain size distribution'], filenum=params['SimFile'], SimDir=params['SimDir'], load_prev=params['Previous Length distribution'], lower_lim=params['Lower grain size limit'], upper_lim=params['Upper grain size limit']):
    if load_prev != None:
        L_array = np.loadtxt(str(f'{path_here}\\{load_prev}'))
        dimensions = L_array.shape
        if (dimensions[0]!=n) and (dimensions[1]!=2*m-2):
            print("WARNING: Loaded length array does not match requested dimensions.")
    elif dist_type=="exponential" or dist_type=="Exponential":
        L = np.random.exponential(L_mean, size=(n,m))
        L_array = np.zeros((n, 2*m-2))
        L_array[:,0] = L[:,0]
        L_array[:,-1] = L[:,-1]
        grain_counter = 1 
        for i in range(1,2*m-3,2):
            L_array[:,i:i+2] = np.column_stack((L[:,grain_counter],L[:,grain_counter]))
            grain_counter+=1
            # print(str(f'DEBUG: i={i}'))
        #save length distrubution
        L_array = np.clip(L_array, lower_lim, upper_lim)
        np.savetxt(str(f'{path_here}\\{SimDir}\\exp_length_dist_{filenum}.txt'), L_array)
    elif dist_type == "constant" or dist_type == "Constant":
        L_array = np.full((n, 2*m-2), L_mean)
    print(str(f'Length array = {L_array}'))
    return L_array

def resistance_array(L_array, ND=params['ND'], ND_array = None, type=params['Doping type'], n=params['Rows of grains (n)'], m=params['Columns of grains (m)'], L=params['Length'], size_type=params['Grain size distribution']):
    if (type == "constant" or type == "Constant") and (size_type == "constant" or size_type == "Constant"): #only calculate resistance once and fill an array (otherwise calculate each resistance individually)
        R_full_grain = bulk_resistance(W=L)
        
        R_half_grain = bulk_resistance(W=L, L=L/2)
        
        R_half_array = np.full((n,2*m-4),R_half_grain)
        
        R_full_array = np.full((n,1),R_full_grain)
        
        R=np.concatenate((R_full_array,R_half_array,R_full_array),axis=1)
        # # DEBUG:
        # print(str(R))
        print(str(f'R_full_grain = {R_full_grain}'))
        print(str(f'R_half_grain = {R_half_grain}'))
        print(str(f'R_full_array = {R_full_array}'))
        print(str(f'R_half_array = {R_half_array}'))
    elif type == 'alternating' or type == 'Alternating' or size_type != "constant" or size_type != "Constant":
        R = np.zeros((n, 2*m-2))
        for k in range(1, 2*m-3): #columns ignoring first and last
            for i in range(n): #rows
                print(f'R={R[i][k]}')
                R[i][k] = bulk_resistance(W=L_array[i][k], ND=ND_array[i][k], L=L_array[i][k]/2)
                print(f'R={R[i][k]},W={L_array[i][k]},L={L_array[i][k]/2}')
        for i in range(n):
            R[i][0] = bulk_resistance(W=L_array[i][0], ND=ND_array[i][0], L=L_array[i][0])
            R[i][2*m-3] = bulk_resistance(W=L_array[i][2*m-3], ND=ND_array[i][2*m-3], L=L_array[i][2*m-3])
        print(str(f'Resistance array = {R}'))
    return R

def ioSetup(in_loc=params['Input locations'], V=params['Voltage inputs'], gnd_loc=params['Ground locations'], VA=params['VA'], n=params['Rows of grains (n)']):
    # Convert locations into lists (values can be tuples (containing lists (containing constants or ranges) or constants), edge or 0)
    # if type(in_loc['left'])==tuple:
    #     if type(in_loc["left"][0])==list:
    #         for i in range(len(in_loc['left'][0])):
    #             if type(in_loc['left'][0][i]) == range:
    #                 in_left += list(in_loc['left'][0][i])
    #             else: #assuming a number
    #                 in_left += [in_loc['left'][0][i]]
    #     else: #one position only
    #         in_left = in_loc['left'][0]
    # elif in_loc['left'] == "edge":
    #     in_left = list(range(n))
    # elif in_loc['left'] == 0:
    #     in_left = []
    #     print("WARNING: No voltage applied to the left edge")
    
    # Attempt 2
    # Generate node dictionary
    V_node_dict = {}
    for i in V:
        V_node_dict[i] = str(f'N{i}')

    in_loc_left = in_loc["left"]
    if type(in_loc_left) == tuple:
        #Generate list of positions
        in_left = []
        V_name_left = []
        V_val_left = []
        V_node_left = []
        resistor_skip = []
        for i in range(len(in_loc_left[0])):
            if type(in_loc_left[0][i]) == range:
                in_left += list(in_loc_left[0][i])
                V_name_left+=[in_loc_left[1][i]]*len(in_loc_left[0][i])
                V_node_left+=[V_node_dict[in_loc_left[1][i]]]*len(in_loc_left[0][i])
                for node in range(in_loc_left[0][i].start, in_loc_left[0][i].stop-1): #decide which resistors to skip (if length is 1, no resistors will be skipped)
                    resistor_skip.append(node)
            else: 
                in_left += [in_loc_left[0][i]]
                V_name_left.append(in_loc_left[1][i])
                V_node_left.append(V_node_dict[in_loc_left[1][i]])
        for name in V_name_left:
            if name in V_name_left:
                V_val_left.append(V[name])
            else:
                print(str(f'ERROR: {name} not defined in voltage inputs.'))
        print(str(f'Voltage input location (left) = {in_left}\nCorresponding source names = {V_name_left}\nCorresponding source values = {V_val_left}\nCorresponding node names = {V_node_left}\nResistors skipped (top node) = {resistor_skip}'))
        if len(in_left) != len(V_name_left) and len(in_left) != len(V_val_left):
            print("ERROR: Voltage inputs don't match up with locations")

        V_left = [in_left, V_node_left, resistor_skip]
        #Generate corresponding list

    elif in_loc_left == "edge" or in_loc_left == "Edge":
        in_left = list(range(n))
        V_applied = next(iter(V_node_dict.values()))
        V_node_left = [V_applied]*n
        resistor_skip = list(range(n-1))
        print(str(f'Voltage input location (left) = {in_left}\nCorresponding node names = {V_node_left}\nResistors skipped (top node) = {resistor_skip}'))
        V_left = [in_left, V_node_left, resistor_skip]
    # Allocate names to the inputs
    # Decide which to sweep
    # Convert list into nodes (more complicated for top and bottom)
    # Check what resistors can be removed (if one voltage applied over a range)
    # Check that there's no overlap between inputs and outputs and that they're within the size of the sim
    # Check that there's the right number of voltages and corresponding values
        
    out_loc_right = gnd_loc['right']
    if type(out_loc_right) == list:
        out_right = []
        resistor_skip_right = []
        for i in range(len(out_loc_right)):
            if type(out_loc_right[i]) == range:
                out_right += list(out_loc_right[i])
                for node in range(out_loc_right[i].start, out_loc_right[i].stop-1): #decide which resistors to skip (if length is 1, no resistors will be skipped)
                    resistor_skip_right.append(node)
            else: 
                out_right += [out_loc_right[i]]
    elif out_loc_right == "edge" or out_loc_right == "Edge":
        out_right = list(range(n))
        resistor_skip_right = list(range(n-1))
    
    gnd_right = [out_right, resistor_skip_right]
    print(str(f'Voltage output location (right) = {out_right}\nResistors skipped (top node) = {resistor_skip_right}'))

    return V_left, gnd_right

def NetlistWriter(ND, R, L, V_left, gnd_right, title=params['LUT_filenum'], V_in=params['Voltage inputs'], SimDir=params['SimDir'], SimFile=params['SimFile'], n=params['Rows of grains (n)'], m=params['Columns of grains (m)'], Vrange=params['VA'], params_length=len(reduced_params), filenum=params['LUT_filenum'], res=params['Voltage sweep step']):
    path = path_here+'\\'+SimDir
    if not os.path.exists(path):
        os.mkdir(path)
            
    ff = open(path+'\\'+SimFile+'.net', 'w')
    ff.write(f'* {title}\n\n')
    for i in V_in:
        ff.write(str(f'{i} N{i} 0 {V_in[i]}\n'))
    # ff.write('V1 x00y00 0 1\n')

    grain_boundary_seq = np.arange(1,3*m-4, 3) # grain boundary row index sequence 

    for i in range(n):
        grain_counter = 0 # so you can index through resistor array
        for k in range(3*m-3):
            if k == 3*m-4 and i in gnd_right[0]: # if it's the bottom right grounded resistor
                ff.write(str(f'Rx{k}y{i} x{k:02.0f}y{i:02.0f} 0 '))
            elif k == 0 and (i in V_left[0]):#LEFT input location matches with list of inputs
                V_indx = V_left[0].index(i) #assumes NO REPEATS/overlap
                ff.write(str(f'Rx{k}y{i} {V_left[1][V_indx]} x{k+1:02.0f}y{i:02.0f} '))
            else:  
                ff.write(str(f'Rx{k}y{i} x{k:02.0f}y{i:02.0f} x{k+1:02.0f}y{i:02.0f} '))
            # print(str(f'grain_counter={grain_counter}'))
            if k in grain_boundary_seq:
                # do grain boundary things
                ff.write(str(f'R=table(V(x{k:02.0f}y{i:02.0f})-V(x{k+1:02.0f}y{i:02.0f})'))
                print(str(f'i={i},k={k}\n'))
                if not LUTchecker(ND1=ND[i][grain_counter-1], ND2=ND[i][grain_counter], L1=L[i][grain_counter-1], L2=L[i][grain_counter]):
                    LUTwriter(ND1=ND[i][grain_counter-1], ND2=ND[i][grain_counter], L1=L[i][grain_counter-1], L2=L[i][grain_counter])
                with open(path+"\\Lookup Tables\\"+'\\LUT_'+str(filenum)+"_ND_"+str(f'{ND[i][grain_counter-1]:.2e}_{ND[i][grain_counter]:.2e}')+"_L_"+str(f'{min(L[i][grain_counter-1], L[i][grain_counter]):.2e}')+'.txt', 'r') as datafile:
                    print("Reading: "+path+"\\Lookup Tables\\"+'\\LUT_'+str(filenum)+"_ND_"+str(f'{ND[i][grain_counter-1]:.2e}_{ND[i][grain_counter]:.2e}')+"_L_"+str(f'{min(L[i][grain_counter-1], L[i][grain_counter]):.2e}')+'.txt')
                    LUTtable = csv.reader(datafile, delimiter=' ')
                    for _ in range(params_length+1):
                        # only starts writing after parameters and column headers
                        next(LUTtable)
                    for ROWS in LUTtable:
                        ff.write(','+str(ROWS[0])+','+str(ROWS[1]))
                ff.write(')\n')
            else:
                ff.write(str(f'{R[i][grain_counter]}\n'))
                grain_counter+=1
    for k in range(0,3*m-2,3): # every column except last (done separately bc of awkward indexing and need for ground)
        grain_counter = 0 #allow you to iterate through ND as it's a nx2m-2 matrix (which doesn't correspond to the nodes)
        for i in range(n-1): # every row in the column (indexes to penultimate node in the column)
            if (i not in V_left[2] or k!=0) and (i not in gnd_right[1] or k!=3*m-3): #checks whether to skip resistor
                if k == 3*m-3 and (i in gnd_right[0] and i+1 in gnd_right[0]): #this shouldn't happen as the resistor should have been skipped
                    print("WARNING: Please use ranges to specify adjacent ground nodes.")
                    continue
                elif k == 3*m-3 and i in gnd_right[0]:
                    ff.write(str(f'Rx{k}y{i}V 0 x{k:02.0f}y{i+1:02.0f} R=table(0-V(x{k:02.0f}y{i+1:02.0f})'))
                elif k == 3*m-3 and i+1 in gnd_right[0]:
                    ff.write(str(f'Rx{k}y{i}V x{k:02.0f}y{i:02.0f} 0 R=table(V(x{k:02.0f}y{i:02.0f})'))
                elif k==0 and i in V_left[0] and i+1 in V_left[0]: #check whether both nodes are connected to voltage source
                    V_indx1 = V_left[0].index(i)
                    V_indx2 = V_left[0].index(i+1)
                    if V_left[1][V_indx1] != V_left[1][V_indx2]:
                        ff.write(str(f'Rx{k}y{i}V {V_left[1][V_indx1]} {V_left[1][V_indx2]} R=table(V({V_left[1][V_indx1]})-V({V_left[1][V_indx2]})'))
                    else:
                        print("WARNING: Please use ranges to specify adjacent nodes connected to the same source.")
                        continue
                elif k==0 and i in V_left[0]: #check whether top node is connected to voltage source
                    V_indx = V_left[0].index(i)
                    ff.write(str(f'Rx{k}y{i}V {V_left[1][V_indx]} x{k:02.0f}y{i+1:02.0f} R=table(V({V_left[1][V_indx]})-V(x{k:02.0f}y{i+1:02.0f})'))
                elif k==0 and i+1 in V_left[0]: #check whether node is connected to voltage source
                    V_indx = V_left[0].index(i+1)
                    ff.write(str(f'Rx{k}y{i}V x{k:02.0f}y{i:02.0f} {V_left[1][V_indx]} R=table(V(x{k:02.0f}y{i:02.0f})-V({V_left[1][V_indx]})'))
                else:
                    ff.write(str(f'Rx{k}y{i}V x{k:02.0f}y{i:02.0f} x{k:02.0f}y{i+1:02.0f} R=table(V(x{k:02.0f}y{i:02.0f})-V(x{k:02.0f}y{i+1:02.0f})'))
                if not LUTchecker(ND1=ND[i][grain_counter], ND2=ND[i+1][grain_counter], L1=L[i][grain_counter], L2=L[i+1][grain_counter]):
                        LUTwriter(ND1=ND[i][grain_counter], ND2=ND[i+1][grain_counter], L1=L[i][grain_counter], L2=L[i+1][grain_counter])
                with open(path+"\\Lookup Tables\\"+'\\LUT_'+str(filenum)+"_ND_"+str(f'{ND[i][grain_counter]:.2e}_{ND[i+1][grain_counter]:.2e}')+"_L_"+str(f'{min(L[i][grain_counter], L[i+1][grain_counter]):.2e}')+'.txt', 'r') as datafile:
                    print('Reading: '+path+"\\Lookup Tables\\"+'\\LUT_'+str(filenum)+"_ND_"+str(f'{ND[i][grain_counter]:.2e}_{ND[i+1][grain_counter]:.2e}')+"_L_"+str(f'{min(L[i][grain_counter], L[i+1][grain_counter]):.2e}')+'.txt')
                    LUTtable = csv.reader(datafile, delimiter=' ')
                    for _ in range(params_length+1):
                        # only starts writing after parameters and column headers
                        next(LUTtable)
                    for ROWS in LUTtable:
                        ff.write(','+str(ROWS[0])+','+str(ROWS[1]))
                ff.write(')\n')
        grain_counter += 1

    ff.write('\n')
    ff.write(str(f'.dc V1 {Vrange[0]} {Vrange[1]} {res}\n'))
    ff.write(".end")
    ff.close()

def LUTchecker(SimDir=params['SimDir'], filenum=params['LUT_filenum'], ND1=1e16, ND2=1e16, L1=1e-4, L2=1e-4, parameters=reduced_params):
    LUT_path = path_here+'\\'+SimDir+'\\Lookup Tables'
    L = min(L1, L2)
    if os.path.exists(LUT_path+'\\LUT_'+str(filenum)+"_ND_"+str(f'{ND1:.2e}_{ND2:.2e}')+"_L_"+str(f'{L:.2e}')+'.txt'):
        with open(LUT_path+'\\LUT_'+str(filenum)+"_ND_"+str(f'{ND1:.2e}_{ND2:.2e}')+"_L_"+str(f'{L:.2e}')+'.txt', 'r', newline='') as csvfile:
            LUTreader = csv.reader(csvfile)
            LUT_params = {}
            for row in csvfile:
                if row and row[0].startswith('#'):
                    # # print(str(f'row = {row}')) 
                    key, value = row[1:].split(':', 1) 
                    # print(str(f'key={key}; value={value}; value_type = {type(value)}'))
                    # LUT_params[key.strip()] = ast.literal_eval(value.strip())
                    LUT_params[key.strip()] = value.strip()
        
        #Set up string dictionary comparison
        string_params = {}
        for p in parameters:
            if type(parameters[p]) == str:
                string_params[p.strip()] = repr(parameters[p]).strip()
            else:
                string_params[p.strip()] = str(parameters[p]).strip()
            # for row in LUTreader:
            #     if row and row[0].startswith('#'):
            #         print(str(f'row = {row}'))
            #         key, value = row[0][1:].split(':', 1)
            #         print(str(f'key={key}, value={value}, value_type = {type(value)}'))
            #         LUT_params[key.strip()] = ast.literal_eval(value.strip())
        params_match = LUT_params==string_params
        # print(str(f'DEBUG: params_match = {params_match}'))
        # print(str(f'DEBUG: table params from existing file = {LUT_params}'))
    else:
        params_match = False
    # print(str(f'DEBUG: simulation params = {string_params}'))
    if not params_match:
        print(str(f'NEW FILE CREATED: {LUT_path}/LUT_{filenum}_ND_{ND1:.2e}_{ND2:.2e}_L_{L:.2e}.txt'))
    return params_match



def AssembleNetlist():
    ND = doping_array()
    R = resistance_array()
    NetlistWriter(ND=ND, R=R)


def Simulate(SimDir, SimFile):
    os.popen("\"C:\\Program Files\\LTC\\LTspiceXVII\\XVIIx64.exe\" -b -ascii \"" + path_here+'\\'+SimDir+'\\'+SimFile+'.net'"\"")
    time.sleep(5)

path = path_here+'\\'+params['SimDir']
if not os.path.exists(path):
    os.mkdir(path)
if params['Contact type']=="point":
    V_left, gnd_right = ioSetup(in_loc={"left":([0],["V1"]), "top":None, "right":None, "bottom":None}, V={"V1":1}, gnd_loc={"left":0, "top":0, "right":[params['Rows of grains (n)']-1], "bottom":None})
else:
    V_left, gnd_right = ioSetup()
ND = doping_array()
L = length_array()
R = resistance_array(L_array=L, ND_array=ND)
NetlistWriter(ND=ND, R=R, L=L, V_left=V_left, gnd_right=gnd_right, title=title)
Simulate(SimDir=params["SimDir"], SimFile=params["SimFile"])

