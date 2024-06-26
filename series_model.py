from sympy import *
import numpy as np
import matplotlib.pyplot as plt
import csv
import os
import pandas as pd
import ltspice
import time
import math

path_here = os.path.dirname(os.path.realpath(__file__))

params = {"ND": [1e16, 2e16], #can be less than the number of GBs and it will repeat the list, or do [1e16] for constant doping
          "VA": (-0.9, 0.9),
          "NT": 4e11, #default 5e11
          "Area": 1.8e-9, #normal size: 1e-9, 3e-11, 4e-10
          "Length": 1e-4, #normal: 1e-4, 2e-5
          "Number of boundaries": 1,
          "SimDir": 'Mobility tests Term 2',
          "SimFile": 'NT 4e11 ND 1e16 2e16 n 1 0.18x1x1 um',
          "LUT_filenum": 'NT 4e11 ND 1e16 2e16 n 1 0.18x1x1 um'}

def find_widths(T=300, er=11.7, ND1=1e16, ND2=1e16, NT=5e11, VA=0): #finds depletion width
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
    else:
        x1s = float(result[1][0])
        x2s = float(result[1][1])
    # print("x1 = "+str(x1s)+" and x2 = "+str(x2s)+"\n")

    return x1s, x2s

def find_barriers(x1, x2, T=300, ND1=1e16, ND2=1e16, NC=2e19, er=11.7):
    q = 1.6e-19
    e0 = 8.85e-14
    k = 8.6e-5
    kT = k*T
    ee = e0*er

    phi_1 = q*ND1*x1**2/(2*ee) + kT*np.log(NC/ND1) #calculating Ec - Efx
    phi_2 = q*ND2*x2**2/(2*ee) + kT*np.log(NC/ND2)  
    return phi_1, phi_2

def calc_current(phi_1, phi_2, A=1, T=300, AA=120):
    k = 8.6e-5
    kT = k*T
    q = 1.6e-19
    current = A*AA*(T**2)*(np.exp(-phi_1/kT) - np.exp(-phi_2/kT))
    return current

def LUTwriter(SimDir, filenum, LUT_iter, ND1=1e16, ND2=1e16, NT=5e11, Vrange=(-1, 1), er=11.7, T=300, res=100, NC=2e19, A=1, AA=130.8):
    path = path_here+'\\'+SimDir
    if not os.path.exists(path):
        os.mkdir(path)
    LUT_path = path+'\\Lookup Tables'
    if not os.path.exists(LUT_path):
        os.mkdir(LUT_path)
    VA = np.linspace(Vrange[0], Vrange[1], res) #makes array from -1 to 1 with 100 values (res for resolution)
    I = np.zeros_like(VA)
    for k in range(len(VA)):
        print("VA="+str(VA[k])+": ")
        x1, x2 = find_widths(VA=VA[k], ND1=ND1, ND2=ND2, NT=NT, T=T, er=er)
        phi_1, phi_2 = find_barriers(x1, x2, ND1=ND1, ND2=ND2, T=T, er=er, NC=NC)
        print("phi_1="+str(phi_1)+" phi_2="+str(phi_2))
        I[k] = calc_current(phi_1, phi_2, A=A, T=T, AA=AA)
    R = VA/I
    with open(LUT_path+'\\LUT_'+str(filenum)+"_"+str(LUT_iter)+'.txt', 'w', newline='') as csvfile:
        Lwriter = csv.writer(csvfile, delimiter=' ')
        # Lwriter.writerow(['VA', 'R'])
        for k in range(len(VA)):
            Lwriter.writerow([str(VA[k]), str(R[k]), str(I[k])])

def bulk_resistance(ND=[1e16], A=1, L=1, n=3):
    q =  1.6e-19
    ND_len = len(ND)
    print(str(f'No. of doping concs = {ND_len}'))
    R = []
    if ND_len == 1:
        for k in range(n+1):
            # R[k] = L/(q*A*(ND[0]*mu_n))
            mu_mono=1400/math.sqrt(1+(ND[0]/3e16)*350/(ND[0]/3e16+350))
            print(f'mu_mono = {mu_mono}')
            R.insert(k, L/(q*A*(ND[0]*mu_mono)))
    elif ND_len <(n+1):#number of doping concs less than number of grains
        i = 0
        for k in range(n+1):
           if i == ND_len:
                i=0
           # R[k] = L/(q*A*(ND[i]*mu_n))
           mu_mono=1400/math.sqrt(1+(ND[i]/3e16)*350/(ND[i]/3e16+350))
           print(f'mu_mono = {mu_mono}')
           R.insert(k, L/(q*A*(ND[i]*mu_mono)))
           i += 1
    else:
        for k in range(n+1):
            # R[k] = L/(q*A*(ND[k]*mu_n))
            mu_mono=1400/math.sqrt(1+(ND[k]/3e16)*350/(ND[k]/3e16+350))
            print(f'mu_mono = {mu_mono}')
            R.insert(k, L/(q*A*(ND[k]*mu_mono)))
    return R

def NetlistWriter(SimDir, SimFile, filenum, n, R, LUT_quant, Vrange=(-1,1)):
    path = path_here+'\\'+SimDir
    if not os.path.exists(path):
        os.mkdir(path)

    #LUT iteration list
    LUT_num = []
    if LUT_quant < n :
        counter = 1
        for i in range(n):
            if counter<=LUT_quant:
                # LUT_num[i] = counter
                LUT_num.insert(i, counter)
                counter += 1
            else:
                counter = 1
                LUT_num.insert(i, counter)
                # LUT_num[i] = counter
                counter += 1
            
    ff = open(path+'\\'+SimFile+'.net', 'w')
    ff.write('* The title\n\n')
    ff.write('V1 x00 0 1\n')

    boundary_counter = 0

    for k in range(2*n+1):
        if k == 2*n:
            ff.write('R'+str(k)+' x'+str(f'{k:02.0f} 0 ')+str(R[int(k/2)])+'\n')
        elif k%2 == 0:
            ff.write('R'+str(k)+' x'+str(f'{k:02.0f} x{k+1:02.0f} ')+str(R[int(k/2)])+'\n')
        elif LUT_quant == 1:
            ff.write('R'+str(k)+' x'+str(f'{k:02.0f} x{k+1:02.0f} ')+'R=table(V(x'+str(f'{k:02.0f}')+')-V(x'+str(f'{k+1:02.0f}')+')')
            with open(path+"\\Lookup Tables\\"+"LUT_"+str(filenum)+"_"+str(LUT_quant)+".txt", 'r') as datafile:
                LUTtable = csv.reader(datafile, delimiter=' ')
                for ROWS in LUTtable:
                    ff.write(','+str(ROWS[0])+','+str(ROWS[1]))
            ff.write(')\n')
        elif LUT_quant < n:
            ff.write('R'+str(k)+' x'+str(f'{k:02.0f} x{k+1:02.0f} ')+'R=table(V(x'+str(f'{(k):02.0f}')+')-V(x'+str(f'{k+1:02.0f}')+')')
            with open(path+"\\Lookup Tables\\"+"LUT_"+str(filenum)+"_"+str(LUT_num[boundary_counter])+".txt", 'r') as datafile:
                LUTtable = csv.reader(datafile, delimiter=' ')
                for ROWS in LUTtable:
                    ff.write(','+str(ROWS[0])+','+str(ROWS[1]))
            ff.write(')\n')
            boundary_counter += 1
        elif LUT_quant >= n:
            ff.write('R'+str(k)+' x'+str(f'{k:02.0f} x{k+1:02.0f} ')+'R=table(V(x'+str(f'{(k):02.0f}')+')-V(x'+str(f'{k+1:02.0f}')+')')
            with open(path+"\\Lookup Tables\\"+"LUT_"+str(filenum)+"_"+str(boundary_counter+1)+".txt", 'r') as datafile:
                LUTtable = csv.reader(datafile, delimiter=' ')
                for ROWS in LUTtable:
                    ff.write(','+str(ROWS[0])+','+str(ROWS[1]))
            ff.write(')\n')
            boundary_counter += 1


    ff.write('\n')
    ff.write('.dc V1 '+str(Vrange[0])+' '+str(Vrange[1])+' 0.01\n')
    ff.write(".print(I(R1))\n")
    ff.write(".end")
    ff.close()

def AssembleNetlist(SimDir, SimFile, filenum="test", A=1, L=1, mu_n=1500, n=3, ND=[1e16], NT=5e11, Vrange=(-1, 1), er=11.7, T=300, res=100, NC=2e19, AA=130.8):
    ND_len = len(ND)
    if ND_len ==1:
        LUT_iter = 1
        LUTwriter(SimDir, filenum, LUT_iter=LUT_iter, ND1=ND[0], ND2=ND[0], NT=NT, Vrange=Vrange, er=er, T=T, res=res, NC=NC, A=A, AA=AA) 
    elif ND_len < (n+1): #need n+1 because there are n+1 grains (where n is the grain boundaries) 
        for k in range(ND_len-1):
            LUT_iter = k + 1
            LUTwriter(SimDir, filenum, LUT_iter=LUT_iter, ND1=ND[k], ND2=ND[k+1], NT=NT, Vrange=Vrange, er=er, T=T, res=res, NC=NC, A=A, AA=AA)
        #loop to start to generate the final table
        LUT_iter=ND_len
        LUTwriter(SimDir, filenum, LUT_iter=LUT_iter, ND1=ND[ND_len-1], ND2=ND[0], NT=NT, Vrange=Vrange, er=er, T=T, res=res, NC=NC, A=A, AA=AA)
    else:
        for k in range(n):
            LUT_iter = k + 1
            LUTwriter(SimDir, filenum, LUT_iter=LUT_iter, ND1=ND[k], ND2=ND[k+1], NT=NT, Vrange=Vrange, er=er, T=T, res=res, NC=NC, A=A, AA=AA)
    R = bulk_resistance(ND=ND, A=A, L=L, n=n)
    NetlistWriter(SimDir=SimDir, SimFile=SimFile, filenum=filenum, n=n, R=R, LUT_quant=LUT_iter, Vrange=Vrange)


AssembleNetlist(SimDir=params["SimDir"], SimFile=params["SimFile"], filenum=params['LUT_filenum'], n=params['Number of boundaries'], ND=params['ND'], Vrange=params["VA"], NT=params['NT'], A=params['Area'], L=params['Length'])

def Simulate(SimDir, SimFile):
    os.popen("\"C:\\Program Files\\LTC\\LTspiceXVII\\XVIIx64.exe\" -b -ascii \"" + path_here+'\\'+SimDir+'\\'+SimFile+'.net'"\"")
    time.sleep(5)

Simulate(SimDir=params["SimDir"], SimFile=params["SimFile"])

