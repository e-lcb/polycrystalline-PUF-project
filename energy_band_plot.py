from sympy import *
import numpy as np
import matplotlib.pyplot as plt
import csv
import os

path_here = os.path.dirname(os.path.realpath(__file__))

def find_widths(T=300, er=11.7, ND1=1e16, ND2=1e16, NT=5e12, VA=0): #finds depletion width
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
    print(str(result))
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

def calc_bands(x1s, x2s, L=1e-4, ND1=1e16, ND2=1e16, NC=2e19, T=300, er=11.7, Eg = 1.12, res=100, VA=0):
    k = 8.6e-5
    kT = k*T
    q = 1.6e-19
    e0 = 8.85e-14
    ee = e0*er
    # #Voltage divider based on bulk resistances?
    # VA1 = VA * R1/(R1+R2)
    # VA2 = VA - VA1
    #Calculate bands
    # x_val = np.linspace(2*-x1s, 2*x2s, res)
    x_val = np.linspace(-L/2, L/2, res)
    cond_band = np.zeros_like(x_val)
    val_band = np.zeros_like(x_val)
    fermi_lvl = np.zeros_like(x_val)
    fermi_lvl_0 = np.zeros_like(x_val)
    for k in range(len(x_val)):
        if x_val[k] < -x1s:
            cond_band[k]=kT*np.log(NC/ND1) 
        elif -x1s <= x_val[k] <= 0:
            cond_band[k]=kT*np.log(NC/ND1) + q*ND1*(x_val[k]+x1s)**2/(2*ee) 
        elif 0 < x_val[k] <= x2s:
            cond_band[k]=kT*np.log(NC/ND2) + q*ND2*(x_val[k]-x2s)**2/(2*ee) - VA
            # +kT * np.log(ND1/ND2) - VA
        elif x_val[k] > x2s:
            cond_band[k]=kT*np.log(NC/ND2) - VA
            # +kT * np.log(ND1/ND2) - VA
        else:
            print("WARNING: Number not in range")
            x_val[k] = 0
        if x_val[k] == 0:
            print(f'Conduction boundary energy = {cond_band[k]}')
        val_band[k] = cond_band[k] - Eg
        fermi_lvl[k]=cond_band[k] - kT*np.log(NC/ND1)
        # fermi_lvl_0[k]=0
    return cond_band, val_band, fermi_lvl, fermi_lvl_0, x_val

def calc_bulk_res(ND1=1e16, ND2=1e16, A=1, L=1, mu_n=1500, mu_p=450): #needs to be adjusted for p-type (assumes ND>>ni)
    q =  1.6e-19
    R1 = L/(q*A*(ND1*1500))
    R2 = L/(q*A*(ND2*1500))
    return R1, R2

def LUTwriter(path, filenum, ND1=1e16, ND2=1e16, NT=5e12, Vrange=(-1, 1), er=11.7, T=300, res=100, NC=2e19, A=1, AA=120):
    VA = np.linspace(Vrange[0], Vrange[1], res) #makes array from -1 to 1 with 100 values (res for resolution)
    I = np.zeros_like(VA)
    for k in range(len(VA)):
        x1, x2 = find_widths(VA=VA[k], ND1=ND1, ND2=ND2, NT=NT, T=T, er=er)
        phi_1, phi_2 = find_barriers(x1, x2, ND1=ND1, ND2=ND2, T=T, er=er, NC=NC)
        I[k] = calc_current(phi_1, phi_2, A=A, T=T, AA=AA)
    R = VA/I
    with open(path+'LUT_'+str(filenum)+'.txt', 'w', newline='') as csvfile:
        Lwriter = csv.writer(csvfile, delimiter=' ')
        # Lwriter.writerow(['VA', 'R'])
        for k in range(len(VA)):
            Lwriter.writerow([str(VA[k]), str(R[k])])

# LUTwriter("""C:\\Users\\elc20\\OneDrive\\Documents\\University\\Year 4\\Project\\Lookup Tables\\""", 2)

params = {"ND1": 1e16,
          "ND2": 1e16,
          "VA": -1.23,
          "L": 1e-4,
          "NT": 4e11,
          "Voltage split file": "voltage_split_data_doping.csv",
          "Add to file?": False, 
          }

# voltage list [0,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9]

for i in [1e16]:
    # params['VA']=i
    params["ND2"] = i
    x1, x2 = find_widths(T=300, er=11.7, ND1=params["ND1"], ND2=params["ND2"], NT=params["NT"], VA=params["VA"])
    print("X1="+str(x1)+ "\nx2="+str(x2))
    # R1, R2 = calc_bulk_res(ND1=1e16, ND2=1e16, A=1, L=1, mu_n=1500, mu_p=450)
    Ec, Ev, Ef, Ef0, x_values = calc_bands(x1, x2, L=params["L"], ND1=params["ND1"], ND2=params["ND2"], NC=2e19, T=300, Eg = 1.12, res=1000, VA=params["VA"])

    print(f'Left conduction band = {Ec[0]}, x = {x_values[0]}')
    print(f'Right conduction band = {Ec[-1]}, x = {x_values[-1]}')
    index_closest = np.abs(x_values).argmin()
    print(f'Boundary condtion band = {Ec[index_closest]}, x = {x_values[index_closest]}')

    if params['VA'] == 0:
        phi1 = Ec[index_closest]-Ec[0]
        phi2 = Ec[index_closest]-Ec[-1]
        a1 = "N/A"
        a2 = "N/A"
    else:
        phi1 = "see above"
        phi2 = "see above"
        a1 = "tbd"
        a2 = "tbd"
    data_names = [["ND1", "ND2", "VA", "V1", "V2", "Vb", "V1-Vb", "V2-Vb", "a1", "a2", "x1", "x2", "phi1", "phi2"]]
    voltage_split_data = [[params['ND1'], params['ND2'], params['VA'], Ec[0], Ec[-1], Ec[index_closest], Ec[index_closest]-Ec[0], Ec[index_closest]-Ec[-1], a1, a2, x1, x2, phi1, phi2]]

    filepath = path_here + '\\' + params['Voltage split file']

    if not os.path.exists(filepath):
        print(f"The file '{filepath}' does not exist.")
        with open(filepath, 'a', newline='', encoding='utf-8') as file:
            writer = csv.writer(file, delimiter=';')
            writer.writerows(data_names)

    if params['Add to file?']:
        print("Adding to file...")
        with open(filepath, 'a', newline='', encoding='utf-8') as file:
            writer = csv.writer(file, delimiter=';')
            writer.writerows(voltage_split_data)

    plt.plot(x_values*10**4, Ec, 'b-', label="Conduction band")
    plt.plot(x_values*10**4, Ev, 'r-', label="Valence band")
    plt.plot([-x1*10**4, x_values[index_closest]*10**4, x2*10**4], [Ec[0], Ec[index_closest], Ec[-1]], 'k.')
    # plt.plot(x_values*10**4, Ef0, 'g--', label="Fermi level")
    # plt.plot(x_values*10**4, Ec, 'b-', label="Conduction band", x_values*10**4, Ev, 'r-', label="Valence band", x_values*10**4, Ef0, 'b--', label="Fermi level")
    # plt.title('Energy levels of grain boundary at VA=0.5V ND1 1e16 ND2 2e16')
    plt.xlabel('x ($\mu$m)')
    plt.ylabel('Energy (eV)')
    # plt.ylim(-1, 0.8)
    plt.xlim(-0.5, 0.5)
    # plt.axvline(x=-x1*10**4, label="x1", alpha=0.5)
    # plt.axvline(x=x2*10**4, label="x2", alpha=0.5)
    # plt.axvspan(-x1*10**4, x2*10**4, alpha=0.2, color='red', label='Depletion region')
    plt.axvspan(-x1*10**4, x2*10**4, alpha=0.2, color='red', label=""+str(f'$x_{1}$={-x1*10**4:.3f} $\mu$m and\n $x_{2}$={x2*10**4:.3f} $\mu$m'))
    # plt.plot(x_values, Ec, 'b-')
    plt.legend(loc='lower left')
    plt.show()
