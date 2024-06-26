import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit as cf
import ltspice
from sympy import *


params = {'ND1': 1e16,
          'ND2': 2e16,
          'NT': 4e11,
        #   'VA': (0.01, 0.9),
          'res': 100,
          'Sim file': 'C:\\Users\\elc20\\OneDrive\\Documents\\University\\Year 4\\Project\\1D Model\\Mobility tests Term 2\\NT 4e11 ND 1e16 2e16 n 1 0.18x1x1 um.raw',
          'A': 1.8e-9,
          'L': 1e-4,
          'N': 1} #REMEMBER TO EDIT THIS

# def h(x, a, b, c):
#     return a*(x**2) + b*x + c

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

def find_barriers(x1, x2, T=300, ND1=1e16, ND2=2e16, NC=2e19, er=11.7):
    q = 1.6e-19
    e0 = 8.85e-14
    k = 8.6e-5
    kT = k*T
    ee = e0*er

    phi_1 = q*ND1*x1**2/(2*ee) + kT*np.log(NC/ND1) #calculating Ec - Efx
    phi_2 = q*ND2*x2**2/(2*ee) + kT*np.log(NC/ND2)  
    return phi_1, phi_2


q = 1.6e-19

R = []

filepath = params['Sim file']
LT = ltspice.Ltspice(filepath)
LT.parse()
time = LT.get_time()
I = LT.get_data('I(R1)')
V = LT.get_data('V1')
V_boundary = LT.get_data('V(x01)') - LT.get_data('V(x02)')
print(str(f'Boundary Voltage = {V_boundary}\n'))
R = V/I
R = np.delete(R, 90) #**Removing value where V is closest to zero
n = (params['ND1']+params['ND2'])/2
mu_eff = params['L']*2/(params['A']*q*R*n)
mu_eff_recip = 1/mu_eff

V_grain = np.zeros((params['N']+1,len(I)))
V_grain1 = np.zeros_like(V)
grain_counter = 0
for i in range(params['N']):
    print(str(f'i={i}'))
    V_grain[i]=LT.get_data(str(f'V(x{grain_counter:02.0f})')) - LT.get_data(str(f'V(x{(grain_counter+1):02.0f})'))
    print(str(f'V(x{grain_counter:02.0f})-V(x{(grain_counter+1):02.0f}):{V_grain[i]}'))
    grain_counter += 2
V_grain[params['N']] = LT.get_data(str(f'V(x{grain_counter:02.0f})')) #final voltage is ground
print(str(f'V(x{grain_counter:02.0f}):{V_grain[params["N"]]}'))

V_grain1 = V_grain[i]
print(str(f'V_grain1 = {V_grain1}'))
# print(str(f'Grain voltage = {V_grain}\n'))

# print(str(f'Total voltage = {V_grain[0]+V_grain[1]+V_boundary}'))
# print(str(f'Applied voltage = {V}'))

# #Inlcude both voltage directions
#  def mobility_func(X, A, B):
#     V, phi_1, phi_2 = X
#     return A/V*(np.exp(B*phi_1)-np.exp(B*phi_2))

# #Only works for positive voltage
# def mobility_func(X, A, B):
#     V, phi_1, phi_2 = X
#     return A/V*(np.exp(B*phi_1))

# #Initial model meant to include geometry
# def mobility_func(phi_1, A, B, C):
#     return A + B/np.exp(C*phi_1)

# Include geometry with voltage dependent terms (this one works with netlist edit provided weird values are removed - when VA is almost zero)
# def mobility_func(X, A, B, C, D):
#     V, V_g, x1, x2, phi_1, phi_2 = X
#     return A/V*(B*(np.exp(C*phi_1)-np.exp(C*phi_2))+D*V_g/((1e-4)-x1-x2))

# Include geometry with voltage dependent terms and both grain voltage (this one works with netlist edit provided weird values are removed - when VA is almost zero)
# def mobility_func(X, A, B, C, D):
#     V, V_g1, V_g2, x1, x2, phi_1, phi_2 = X
#     return A/V*(B*(np.exp(C*phi_1)-np.exp(C*phi_2))+D*(V_g1/((1e-4)-x1-x2)+V_g2/((1e-4)-x2)))

# #Remove D as it's removing it with negative voltage
# def mobility_func(X, A, B, C):
#     V, V_g, x, phi_1, phi_2 = X
#     return A/V*(B*(np.exp(C*phi_1)-np.exp(C*phi_2))+2.24*V_g/((1e-4)-2*x))

# Using conduction to put in form in the paper
# def mobility_func(X, A, B, C):
#     V_boundary, x1, x2, phi_1, phi_2 = X
#     mu_eff_recip = A*(params['L']-x1-x2)+B*V_boundary/(np.exp(C*phi_1)-np.exp(C*phi_2))
#     return 1/mu_eff_recip

# Using conduction to put in form in the paper SCALED UP (works mostly)
def mobility_func(X, A, B, C, D):
    V_boundary, x1, x2, phi_1, phi_2 = X
    mu_eff_recip = A*(params['L']-x1)+B*V_boundary/(np.exp(C*phi_1)-np.exp(C*phi_2))+D*(params['L']-x2)
    return 1/mu_eff_recip



phi_1 = np.zeros_like(V)
phi_2 = np.zeros_like(V)
x1 = np.zeros_like(V)
x2 = np.zeros_like(V)

for i in range(len(V)):
    x1[i], x2[i]= find_widths(ND1=params['ND1'], ND2=params['ND2'], VA=V_boundary[i], NT=params['NT']) #not sure this is the right value for VA as it doesn't include the voltage on the grains
    phi_1[i], phi_2[i] = find_barriers(x1[i], x2[i], ND1=params['ND1'], ND2=params['ND2'])

print(str(f'x1 = {x1}'))

for i in range(len(V)):
    print(str(f'i={i}, VA={V[i]}'))
# f, g = cf(mobility_func, (V, phi_1, phi_2), mu_eff, p0=[ 5.96440890e+05, -3.84292774e+01])#where f is parameter values and g is correlation matrix (how well it fits), p0 is starting guess

# f, g = cf(mobility_func, phi_1, mu_eff_recip, p0=[4.12e-4, 5.96440890e+05, -3.84292774e+01])

#Test for this specific case: \\1D Model\\Mobility tests\\NT 5e11 ND 1e16 2e16 n 1.raw
#replace V_grain[0] with V_grain_removed in cf
x1=np.delete(x1, 90)
x2=np.delete(x2, 90)
phi_1=np.delete(phi_1, 90)
phi_2=np.delete(phi_2, 90)
V_grain1_removed=np.delete(V_grain[0], 90)
V_grain2_removed=np.delete(V_grain[1], 90)
V_boundary = np.delete(V_boundary, 90)
V=np.delete(V, 90)

# for i in range(len(V)):
#     print(str(f'i={i}, VA={V[i]}'))

#This one works with the functioning model to include geometry
# f, g = cf(mobility_func, (V, V_grain1_removed, x1, x2, phi_1, phi_2), mu_eff, p0=[ 4.88813851e-02, 1.41007522e+05, -3.87955746e+01, 1.51903309e+00])

#This one works with the functioning model to include geometry and BOTH grain voltages
# f, g = cf(mobility_func, (V, V_grain1_removed, V_grain2_removed, x1, x2, phi_1, phi_2), mu_eff, p0=[ 1.4e-02, 1.08e+07, -3.87956e+01, 3.36e+00])

#This one is meant to have the nicer form but doesn't work yet
# f, g = cf(mobility_func, (V_boundary, x1, x2, phi_1, phi_2), mu_eff, p0=[ 7.14, 1.48e-6, -3.87955746e+01])

#This one is meant to have the nicer form but doesn't work yet
f, g = cf(mobility_func, (V_boundary, x1, x2, phi_1, phi_2), mu_eff, p0=[ 6.185, 1.02e-6, -38.76, 3.457])

# p0=[ 10.7, 2.2e-6, -3.87955746e+01, 5.4]
# Constant doping initial guesses
# p0=[ 4.12, 7.41e-7, -38.76, 4.12]
# 1e16 2e16 doping guesses
p0=[ 6.185, 1.11e-6, -38.76, 3.457]

#initial guesses for updated model: p0=[0.0625,1.08e7,-3.84292774e+01,2.24]

# #Using previously found terms:
# f = [ 3.46552570e-02, 1.71646631e+07, -4.22694969e+01, 2.30660458e+00]

# p0=[675000, -38.6]

# def find_widths(T=300, er=11.7, ND1=1e16, ND2=1e16, NT=5e11, VA=0)
# def find_barriers(x1, x2, T=300, ND1=1e16, ND2=1e16, NC=2e19, er=11.7):
print(str(f))

plt.figure(1)
# plt.plot(V, mu_eff_recip, 'b-', V, 1/mobility_func(X=(V, phi_1, phi_2), A=f[0], B=f[1]), 'r--')
# plt.plot(V, mu_eff_recip, 'b-', V, mobility_func(phi_1, A=f[0], B=f[1], C=f[2]), 'r--')
# This one works with geometry model: 
# plt.plot(V, mu_eff_recip, 'b-', V, 1/mobility_func(X=(V, V_grain1_removed, V_grain2_removed, x1, x2, phi_1, phi_2), A=f[0], B=f[1], C=f[2], D=f[3]), 'r--')
# This one is in paper form but doesn't really work (also take reciprocal for this one)
plt.plot(V, mu_eff_recip, 'b-', label="Simulation")
# plt.plot(V, 1/mobility_func(X=(V_boundary, x1, x2, phi_1, phi_2), A=f[0], B=f[1], C=f[2]), 'r--', label="Model")
plt.plot(V, 1/mobility_func(X=(V_boundary, x1, x2, phi_1, phi_2), A=f[0], B=f[1], C=f[2], D=f[3]), 'r--', label="Model")
plt.plot(V, 1/mobility_func(X=(V_boundary, x1, x2, phi_1, phi_2), A=p0[0], B=p0[1], C=p0[2], D=p0[3]), 'g--', label="Guess")
plt.xlabel('Applied voltage (V)')
plt.ylabel('1/Effective mobility (Vs/cm$^{2}$)')
plt.legend()
plt.figure(2)
# plt.plot(V, mu_eff, 'b-', V, mobility_func(X=(V, phi_1, phi_2), A=f[0], B=f[1]), 'r--')
# plt.plot(V, mu_eff, 'b-', V, mobility_func(phi_1, A=f[0], B=f[1], C=f[2]), 'r--')
# plt.plot(V, mu_eff, 'b-', V, mobility_func(X=(V, V_grain1_removed, V_grain2_removed, x1, x2, phi_1, phi_2), A=f[0], B=f[1], C=f[2], D=f[3]), 'r--')
plt.plot(V, mu_eff, 'b-', label="Simulation")
# plt.plot(V, mobility_func(X=(V_boundary, x1, x2, phi_1, phi_2), A=f[0], B=f[1], C=f[2]), 'r--', label="Model")
plt.plot(V, mobility_func(X=(V_boundary, x1, x2, phi_1, phi_2), A=f[0], B=f[1], C=f[2], D=f[3]), 'r--', label="Model")
plt.plot(V, mobility_func(X=(V_boundary, x1, x2, phi_1, phi_2), A=p0[0], B=p0[1], C=p0[2], D=p0[3]), 'g--', label="Guess")
plt.xlabel('Applied voltage (V)', fontsize=16)
plt.ylabel('Effective mobility (cm$^{2}$V$^{-1}$s$^{-1}$)', fontsize=16)
plt.legend(fontsize=14)

plt.figure(3)
# plt.semilogy(V, mu_eff, 'b-', V, mobility_func(X=(V, phi_1, phi_2), A=f[0], B=f[1]), 'r--')
# plt.semilogy(V, mu_eff, 'b-', V, mobility_func(phi_1, A=f[0], B=f[1], C=f[2]), 'r--')
# plt.semilogy(V, mu_eff, 'b-', V, mobility_func(X=(V, V_grain1_removed, V_grain2_removed, x1, x2, phi_1, phi_2), A=f[0], B=f[1], C=f[2], D=f[3]), 'r--')
plt.semilogy(V, mu_eff, 'b-', label="Simulation")
# plt.semilogy(V, mobility_func(X=(V_boundary, x1, x2, phi_1, phi_2), A=f[0], B=f[1], C=f[2]), 'r--', label="Model")
plt.semilogy(V, mobility_func(X=(V_boundary, x1, x2, phi_1, phi_2), A=f[0], B=f[1], C=f[2], D=f[3]), 'r--', label="Best fit")
plt.semilogy(V, mobility_func(X=(V_boundary, x1, x2, phi_1, phi_2), A=p0[0], B=p0[1], C=p0[2], D=p0[3]), 'g--', label="Expected", alpha=0.8)
# plt.semilogy(V, mobility_func(X=(V_boundary, x1, x2, phi_1, phi_2), A=p0[0], B=p0[1], C=p0[2], D=p0[3]), 'g--', label="Model with expected parameters")

plt.xlabel('Applied voltage (V)', fontsize=16)
plt.ylabel('log[Effective mobility (cm$^{2}$V$^{-1}$s$^{-1}$)]', fontsize=16)
plt.legend(fontsize=14)
plt.show()
