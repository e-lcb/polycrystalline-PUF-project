from sympy import *
import numpy as np
import matplotlib.pyplot as plt
import ltspice

raw_filepath_prefix = 'C:\\Users\\elc20\\OneDrive\\Documents\\University\\Year 4\\Project\\1D Model\\Initial report diagrams\\'

params = {'ND':[1e16],
          '.raw file':raw_filepath_prefix + 'NT 5e11 ND 1e16 n 5 0.18x1x1 um.raw',
          'n':5, 
          'VA':4, #it will find the closest value in the simulation
          'T':300,
          'er':11.7,
          'Eg':1.12,
          'L':1e-4, #1e-4
          'NC':2e19, 
          'NT':5e11,
          'res':5000} #needs to be divisible by n

#note: Change .raw file to simulation with suitable conditions to work out correct voltage split (unless no bias is applied)

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
    
    return x1s, x2s

def find_barriers(x1, x2, T=300, ND1=1e16, ND2=1e16, NC=2e19, er=11.7):
    q = 1.6e-19
    e0 = 8.85e-14
    k = 8.6e-5
    kT = k*T
    ee = e0*er

    phi_1 = q*ND1*x1**2/(2*ee) + kT*np.log(NC/ND1) #calculating Ec - Efx
    phi_2 = q*ND2*x2**2/(2*ee) + kT*np.log(NC/ND2) + kT * np.log(ND1/ND2) 
    return phi_1, phi_2

def calc_current(phi_1, phi_2, A=1, T=300, AA=120):
    k = 8.6e-5
    kT = k*T
    current = A*AA*(T**2)*(np.exp(-phi_1/kT) - np.exp(-phi_2/kT))
    return current

def calc_bands(x1s, x2s, x_val, init_val, res, fermi_lvl_diff, ND1=1e16, ND2=1e16, NC=2e19, T=300, er=11.7, Eg = 1.12, VA=0.5, L=1e-4, n=2):
    k = 8.6e-5
    kT = k*T
    q = 1.6e-19
    e0 = 8.85e-14
    ee = e0*er
    #Calculate bands
    # print('res/n='+str(res/n))
    x_val = np.linspace(-L/2, L/2, int(res/n))
    cond_band = np.zeros_like(x_val)
    val_band = np.zeros_like(x_val) 
    fermi_lvl = np.zeros_like(x_val)
    fermi_lvl0 = np.zeros_like(x_val)
    for k in range(len(x_val)):
        if x_val[k] < -x1s:
            cond_band[k]=kT*np.log(NC/ND1) + init_val
        elif -x1s <= x_val[k] <= 0:
            cond_band[k]=kT*np.log(NC/ND1) + q*ND1*(x_val[k]+x1s)**2/(2*ee) + init_val
        elif 0 < x_val[k] <= x2s:
            cond_band[k]=kT*np.log(NC/ND2) + q*ND2*(x_val[k]-x2s)**2/(2*ee) - VA + init_val
            # +kT * np.log(ND1/ND2) - VA
        elif x_val[k] > x2s:
            cond_band[k]=kT*np.log(NC/ND2) - VA + init_val
            # +kT * np.log(ND1/ND2) - VA
        else:
            print("WARNING: Number not in range")
            x_val[k] = 0
        val_band[k] = cond_band[k] - Eg
        # fermi_lvl[k] = cond_band[k] - kT*np.log(NC/ND1) 
        fermi_lvl[k] = cond_band[k] - fermi_lvl_diff
    init_val = -VA + init_val #so the next grain aligns as the voltage isn't included at the start
    return cond_band, val_band, fermi_lvl, fermi_lvl0, init_val


# LUTwriter("""C:\\Users\\elc20\\OneDrive\\Documents\\University\\Year 4\\Project\\Lookup Tables\\""", 2)

# x1, x2 = find_widths(T=300, er=11.7, ND1=1e16, ND2=2e16, NT=5e11, VA=0.5)
# print("X1="+str(x1)+ "\nx2="+str(x2))
# Ec, Ev, Ef, x_values = calc_bands(x1, x2, ND1=1e16, ND2=2e16, NC=2e19, T=300, Eg = 1.12, res=1000, VA=0.5)

# plt.plot(x_values, Ec, 'b-', x_values, Ev, 'r-', x_values, Ef, 'b--')
# plt.title('Energy levels of grain boundary at VA=0.5V ND1 1e16 ND2 2e16')
# plt.xlabel('x')
# plt.ylabel('Energy (eV)')
# # plt.plot(x_values, Ec, 'b-')
# plt.show()

def LTspice_parser(raw_file=params['.raw file'], n=params['n'], VA=params['VA']):
    filepath = raw_file
    LT = ltspice.Ltspice(filepath)
    LT.parse()
    V_source = LT.get_data('V1')
    V = np.zeros((n, len(V_source)))
    # V = []
    
    indx = np.searchsorted(V_source, VA) #the voltage may not exactly match the desired applied voltage so have to find the closest one
    
    if abs(V_source[indx]-VA) < abs(V_source[indx-1]-VA):
        indx = indx
    else:
        indx = indx - 1
    VA_actual = V_source[indx] #this is the closest value from the simulation
    print('simulation VA='+str(VA_actual)+'\n')

    for i in range(n): #NOT ANYMORE: do the final value separately because the last node is always a 0
        #V is actual voltage across the boundary at each boundary (it's an n x (number of voltage points) array)
        # V[i] = LT.get_data(str(f'V(x{2*i+1:02.0f})-V({2*i+2:02.0f})'))
        V[i] = (LT.get_data(str(f'V(x{2*i+1:02.0f})'))-LT.get_data(str(f'V(x{2*i+2:02.0f})')))
    print(str(f'Voltage across each boundary at all different voltages:\n {V}'))
    # print(str(LT.get_data(f'V(x01)')-LT.get_data(f'V(x02)')))
    return V, indx, VA_actual

def band_diagram_plot(V, indx, VA_actual, ND=[1e16, 2e16], n=2, T=300, er=11.7, L=1e-4, NC=2e19, Eg=1.12, NT=5e11, res=100):

    k = 8.6e-5
    kT = k*T
    # res=len(V[0])

    x_depletion = np.zeros(2*n)
    ND_all = np.zeros(n+1)
    if len(ND) == 1:
        for i in range(n+1):
            ND_all[i] = ND[0]
    elif len(ND) < n+1: #because n is the number of boundaries not grains
        counter = 0
        for i in range(n+1):
            print(str(f'counter1={counter}'))
            if counter > len(ND)-1:
                counter = 0
                print(str(counter))
            ND_all[i] = ND[counter]
            print(str(f'i={i}, counter={counter}\n'))
            counter += 1
    else:
        for i in range(n+1):
            ND_all[i] = ND[i]

    #calculate bands
    x_val = np.linspace(0, n*L, res) #ignore half of first and last grain
    # if res % 2 != 0:
    #     #delete last value of list BUT does resolution from the simulation actually matter??
    #     print('test')
    #     np.delete(x_val, len(x_val)-1)
    cond_band_full = []
    val_band_full = [] 
    fermi_lvl_full = []
    fermi_lvl_full0 = []
    init_val = 0
    fermi_lvl_diff = kT*np.log(NC/ND_all[0])
    
    for i in range(n):
        x_depletion[2*i], x_depletion[2*i+1] = find_widths(ND1=ND_all[i], ND2=ND_all[i+1], VA=V[i][indx], T=T, er=er, NT=NT)
        if (x_depletion[2*i] >= L/2) or (x_depletion[2*i+1] >= L/2):
            print("WARNING: depletion region too wide")
            print('x1=' + str(x_depletion[2*i])+ ', x2=' + str(x_depletion[2*i+1]))
        cond_band_partial, val_band_partial, fermi_lvl_partial, fermi_lvl_partial0, init_val = calc_bands(x_depletion[2*i], x_depletion[2*i+1], x_val, init_val, res=res, fermi_lvl_diff=fermi_lvl_diff, ND1=ND_all[i], ND2=ND_all[i+1], NC=NC, T=T, er=er, Eg=Eg, VA=V[i][indx], L=L, n=n)
        cond_band_full.extend(cond_band_partial)
        val_band_full.extend(val_band_partial)
        fermi_lvl_full.extend(fermi_lvl_partial)
        fermi_lvl_full0.extend(fermi_lvl_partial0)


    print('x value length='+str(len(x_val))+', conduction band length='+str(len(cond_band_full))+'\n')
    print(f'depletion widths: {x_depletion}\n')

    # plt.plot(x_val, cond_band_full, 'b-', x_val, val_band_full, 'r-', x_val, fermi_lvl_full0, 'b--') #for sus fermi level, plot fermi_lvl_full
    plt.plot(x_val*10**4, cond_band_full, 'b-', label="Conduction band")
    plt.plot(x_val*10**4, val_band_full, 'r-', label="Valence band")
    # plt.plot(x_val*10**4, fermi_lvl_full0, 'g--', label="Fermi level")
    # plt.title(str(f'Energy levels of grain boundary at VA={VA_actual:.2f} ND1 {params["ND"][0]:.1e} ND2 {ND_all[1]:.1e}'))
    plt.xlabel('x ($\mu$m)')
    plt.ylabel('Energy (eV)')
    plt.legend()
    # plt.plot(x_values, Ec, 'b-')
    plt.show()


V, indx, VA_actual = LTspice_parser(raw_file=params['.raw file'], n=params['n'], VA=params['VA'])
if params['VA'] == 0:
    V = np.zeros_like(V)
    VA_actual=0

band_diagram_plot(V=V, indx=indx, VA_actual=VA_actual, ND=params['ND'], n=params['n'], er=params['er'], L=params['L'], NC=params['NC'], Eg=params['Eg'], NT=params['NT'], res=params['res'])


    



        

