import numpy as np
import ltspice
import matplotlib.pyplot as plt


filepath_res_only = 'C:\\Users\\elc20\\OneDrive\\Documents\\University\\Year 4\\Project\\simple_model\\Initial report diagrams\\ND1 1e16 ND2 1e16 var res only A 4e-10 L 2e-5.raw'
filepath_two_grains = 'C:\\Users\\elc20\\OneDrive\\Documents\\University\\Year 4\\Project\\simple_model\\Initial report diagrams\\ND1 1e16 ND2 1e16 A 4e-10 L 2e-5.raw'
filepath_two_grains_poisson_extremes = 'C:\\Users\\elc20\\OneDrive\\Documents\\University\\Year 4\\Project\\1D Model\\Initial report diagrams\\NT 5e11 ND 8.25e15 1.19e16 n 1.raw'
filepath_two_grains_poisson_extremes_05x1x1um = 'C:\\Users\\elc20\\OneDrive\\Documents\\University\\Year 4\\Project\\1D Model\\Initial report diagrams\\NT 5e11 ND 9.768e15, 1.0234e16 n 1 1x1x0.5 um.raw'
filepath_two_grains_equal_doping_05x1x1um = 'C:\\Users\\elc20\\OneDrive\\Documents\\University\\Year 4\\Project\\1D Model\\Initial report diagrams\\NT 5e11 ND 1e16 n 1 1x1x0.5 um.raw'
filepath_two_grains_poisson_extremes_018x1x1um = 'C:\\Users\\elc20\\OneDrive\\Documents\\University\\Year 4\\Project\\1D Model\\Initial report diagrams\\NT 5e11 ND 9.5389e15, 1.0467e16 n 1 0.18x1x1 um.raw'
filepath_two_grains_equal_doping_018x1x1um = 'C:\\Users\\elc20\\OneDrive\\Documents\\University\\Year 4\\Project\\1D Model\\Initial report diagrams\\NT 5e11 ND 1e16 n 1 0.18x1x1 um.raw'
filepath_two_grains_1e16_2e16_018x1x1um = 'C:\\Users\\elc20\\OneDrive\\Documents\\University\\Year 4\\Project\\1D Model\\Initial report diagrams\\NT 5e11 ND 1e16 2e16 n 1 0.18x1x1 um.raw'

filepath_three_grains_poisson_extremes_018x1x1um = 'C:\\Users\\elc20\\OneDrive\\Documents\\University\\Year 4\\Project\\1D Model\\Initial report diagrams\\NT 5e11 ND 9.5389e15, 1.0467e16 n 2 0.18x1x1 um.raw'
filepath_three_grains_equal_doping_018x1x1um = 'C:\\Users\\elc20\\OneDrive\\Documents\\University\\Year 4\\Project\\1D Model\\Initial report diagrams\\NT 5e11 ND 1e16 n 2 0.18x1x1 um.raw'
filepath_six_grains_poisson_extremes_018x1x1um = 'C:\\Users\\elc20\\OneDrive\\Documents\\University\\Year 4\\Project\\1D Model\\Initial report diagrams\\NT 5e11 ND 9.5389e15, 1.0467e16 n 5 0.18x1x1 um.raw'
filepath_six_grains_equal_doping_018x1x1um = 'C:\\Users\\elc20\\OneDrive\\Documents\\University\\Year 4\\Project\\1D Model\\Initial report diagrams\\NT 5e11 ND 1e16 n 5 0.18x1x1 um.raw'

## Variable resistor only - Equal doping
LT_var_res = ltspice.Ltspice(filepath_res_only)
LT_var_res.parse()
I_var_res = LT_var_res.get_data('I(R1)')
V_boundary_var_res = LT_var_res.get_data('V(N1)') #should be the same as the applied voltage
V_applied_var_res = LT_var_res.get_data('V1')
## 2 Grains - Equal Doping
LT_2_grains = ltspice.Ltspice(filepath_two_grains)
LT_2_grains.parse()
I_2_grains = LT_2_grains.get_data('I(R1)')
# # Debug
# print(str(f"V(N2)={LT_2_grains.get_data('V(N2)')}"))
# print(str(f"V(N3)={LT_2_grains.get_data('V(N3)')}"))
# print(str(f"I(R1)={LT_2_grains.get_data('I(R1)')}"))
# print(str(f"V1={LT_var_res.get_data('V1')}"))
# print(str(f"I(R1)={LT_var_res.get_data('I(R1)')}"))
V_boundary_2_grains = LT_2_grains.get_data('V(N2)')-LT_2_grains.get_data('V(N3)')
V_applied_2_grains = LT_2_grains.get_data('V1')
## 2 Grains - Poisson extreme doping 
LT_2g_poisson_x = ltspice.Ltspice(filepath_two_grains_poisson_extremes)
LT_2g_poisson_x.parse()
I_2g_poisson_x = LT_2g_poisson_x.get_data('I(R1)')
V_applied_2g_poisson_x = LT_2g_poisson_x.get_data('V1')
V_boundary_2g_poisson_x = LT_2g_poisson_x.get_data('V(x01)') - LT_2g_poisson_x.get_data('V(x02)')
## 2 Grains - Poisson extreme doping 0.5 by 1 by 1 um
LT_2g_poisson_x_05x1x1um = ltspice.Ltspice(filepath_two_grains_poisson_extremes_05x1x1um)
LT_2g_poisson_x_05x1x1um.parse()
I_2g_poisson_x_05x1x1um = LT_2g_poisson_x_05x1x1um.get_data('I(R1)')
V_applied_2g_poisson_x_05x1x1um = LT_2g_poisson_x_05x1x1um.get_data('V1')
V_boundary_2g_poisson_x_05x1x1um = LT_2g_poisson_x_05x1x1um.get_data('V(x01)') - LT_2g_poisson_x_05x1x1um.get_data('V(x02)')
## 2 Grains - Equal doping 0.5 by 1 by 1 um
LT_2g_equal_doping_05x1x1um = ltspice.Ltspice(filepath_two_grains_equal_doping_05x1x1um)
LT_2g_equal_doping_05x1x1um.parse()
I_2g_equal_doping_05x1x1um = LT_2g_equal_doping_05x1x1um.get_data('I(R1)')
V_applied_2g_equal_doping_05x1x1um = LT_2g_equal_doping_05x1x1um.get_data('V1')
V_boundary_2g_equal_doping_05x1x1um = LT_2g_equal_doping_05x1x1um.get_data('V(x01)') - LT_2g_equal_doping_05x1x1um.get_data('V(x02)')
## 2 Grains - Poisson extreme doping 0.18 by 1 by 1 um
LT_2g_poisson_x_018x1x1um = ltspice.Ltspice(filepath_two_grains_poisson_extremes_018x1x1um)
LT_2g_poisson_x_018x1x1um.parse()
I_2g_poisson_x_018x1x1um = LT_2g_poisson_x_018x1x1um.get_data('I(R1)')
V_applied_2g_poisson_x_018x1x1um = LT_2g_poisson_x_018x1x1um.get_data('V1')
V_boundary_2g_poisson_x_018x1x1um = LT_2g_poisson_x_018x1x1um.get_data('V(x01)') - LT_2g_poisson_x_018x1x1um.get_data('V(x02)')
## 2 Grains - Equal doping 0.18 by 1 by 1 um
LT_2g_equal_doping_018x1x1um = ltspice.Ltspice(filepath_two_grains_equal_doping_018x1x1um)
LT_2g_equal_doping_018x1x1um.parse()
I_2g_equal_doping_018x1x1um = LT_2g_equal_doping_018x1x1um.get_data('I(R1)')
V_applied_2g_equal_doping_018x1x1um = LT_2g_equal_doping_018x1x1um.get_data('V1')
V_boundary_2g_equal_doping_018x1x1um = LT_2g_equal_doping_018x1x1um.get_data('V(x01)') - LT_2g_equal_doping_018x1x1um.get_data('V(x02)')
## 2 Grains - 1e16 2e16 0.18 by 1 by 1 um
LT_2g_1e16_2e16_018x1x1um = ltspice.Ltspice(filepath_two_grains_1e16_2e16_018x1x1um)
LT_2g_1e16_2e16_018x1x1um.parse()
I_2g_1e16_2e16_018x1x1um = LT_2g_1e16_2e16_018x1x1um.get_data('I(R1)')
V_applied_2g_1e16_2e16_018x1x1um = LT_2g_1e16_2e16_018x1x1um.get_data('V1')
V_boundary_2g_1e16_2e16_018x1x1um = LT_2g_1e16_2e16_018x1x1um.get_data('V(x01)') - LT_2g_1e16_2e16_018x1x1um.get_data('V(x02)')
## 3 Grains - Poisson extreme doping 0.18 by 1 by 1 um
LT_3g_poisson_x_018x1x1um = ltspice.Ltspice(filepath_three_grains_poisson_extremes_018x1x1um)
LT_3g_poisson_x_018x1x1um.parse()
I_3g_poisson_x_018x1x1um = LT_3g_poisson_x_018x1x1um.get_data('I(R1)')
V_applied_3g_poisson_x_018x1x1um = LT_3g_poisson_x_018x1x1um.get_data('V1')
V_boundary_3g_poisson_x_018x1x1um = LT_3g_poisson_x_018x1x1um.get_data('V(x01)') - LT_3g_poisson_x_018x1x1um.get_data('V(x02)')
## 3 Grains - Equal doping 0.18 by 1 by 1 um
LT_3g_equal_doping_018x1x1um = ltspice.Ltspice(filepath_three_grains_equal_doping_018x1x1um)
LT_3g_equal_doping_018x1x1um.parse()
I_3g_equal_doping_018x1x1um = LT_3g_equal_doping_018x1x1um.get_data('I(R1)')
V_applied_3g_equal_doping_018x1x1um = LT_3g_equal_doping_018x1x1um.get_data('V1')
V_boundary_3g_equal_doping_018x1x1um = LT_3g_equal_doping_018x1x1um.get_data('V(x01)') - LT_3g_equal_doping_018x1x1um.get_data('V(x02)')
## 6 Grains - Poisson extreme doping 0.18 by 1 by 1 um
LT_6g_poisson_x_018x1x1um = ltspice.Ltspice(filepath_six_grains_poisson_extremes_018x1x1um)
LT_6g_poisson_x_018x1x1um.parse()
I_6g_poisson_x_018x1x1um = LT_6g_poisson_x_018x1x1um.get_data('I(R1)')
V_applied_6g_poisson_x_018x1x1um = LT_6g_poisson_x_018x1x1um.get_data('V1')
V_boundary_6g_poisson_x_018x1x1um = LT_6g_poisson_x_018x1x1um.get_data('V(x01)') - LT_6g_poisson_x_018x1x1um.get_data('V(x02)')
## 6 Grains - Equal doping 0.18 by 1 by 1 um
LT_6g_equal_doping_018x1x1um = ltspice.Ltspice(filepath_six_grains_equal_doping_018x1x1um)
LT_6g_equal_doping_018x1x1um.parse()
I_6g_equal_doping_018x1x1um = LT_6g_equal_doping_018x1x1um.get_data('I(R1)')
V_applied_6g_equal_doping_018x1x1um = LT_6g_equal_doping_018x1x1um.get_data('V1')
V_boundary_6g_equal_doping_018x1x1um = LT_6g_equal_doping_018x1x1um.get_data('V(x01)') - LT_6g_equal_doping_018x1x1um.get_data('V(x02)')
## Results plotting
# plt.plot(V_applied_var_res, I_var_res, label="Variable resistor only")
# plt.plot(V_applied_2_grains, I_2_grains*10**9, 'r-', label="Equal doping 0.2x0.2x0.2$\mu$m")
# plt.plot(V_applied_2g_poisson_x, I_2g_poisson_x*10**9, 'b--', label="Poisson extremes 0.2x0.2x0.2$\mu$m")
# plt.plot(V_applied_2g_poisson_x_05x1x1um, I_2g_poisson_x_05x1x1um*10**9, 'm--', label="Poisson extremes 0.5x1x1$\mu$m")
# plt.plot(V_applied_2g_equal_doping_05x1x1um, I_2g_equal_doping_05x1x1um*10**9, 'r-', label="Equal doping 0.5x1x1$\mu$m")
plt.plot(V_applied_2g_equal_doping_018x1x1um, I_2g_equal_doping_018x1x1um*10**6, 'b', label="Equal doping")
plt.plot(V_applied_2g_poisson_x_018x1x1um, I_2g_poisson_x_018x1x1um*10**6, 'r--', label="Poisson extremes")
# plt.plot(V_applied_2g_1e16_2e16_018x1x1um, I_2g_1e16_2e16_018x1x1um*10**6, 'purple', label="1e16 2e16 0.18x1x1$\mu$m")
# plt.plot(V_applied_3g_equal_doping_018x1x1um, I_3g_equal_doping_018x1x1um*10**6, 'g', label="3 grains - Equal doping")
# plt.plot(V_applied_3g_poisson_x_018x1x1um, I_3g_poisson_x_018x1x1um*10**6, 'c--', label="3 grains - Poisson extremes")
# plt.plot(V_applied_6g_equal_doping_018x1x1um, I_6g_equal_doping_018x1x1um*10**6, 'g', label="6 grains - Equal doping") #label="6 grains - Equal doping"
# plt.plot(V_applied_6g_poisson_x_018x1x1um, I_6g_poisson_x_018x1x1um*10**6, 'c--', label="6 grains - Poisson extremes")

plt.xlabel('Applied voltage (V)', fontsize=14)
plt.ylabel('Current ($\mu$A)', fontsize=14)
plt.xlim(-1.5, 1.5)
# plt.ylim(-1.9, 2.3)
plt.legend(loc='lower right', fontsize=11)

# V_applied_2g_equal_doping_018x1x1um = np.delete(V_applied_2g_equal_doping_018x1x1um, 150)
# I_2g_equal_doping_018x1x1um = np.delete(I_2g_equal_doping_018x1x1um, 150)
# V_applied_2g_poisson_x_018x1x1um = np.delete(V_applied_2g_poisson_x_018x1x1um, 150)
# I_2g_poisson_x_018x1x1um = np.delete(I_2g_poisson_x_018x1x1um, 150)

plt.figure()
plt.plot(V_applied_2g_equal_doping_018x1x1um, np.log10(abs(I_2g_equal_doping_018x1x1um)), 'r-', label="2 grains - Equal doping")
plt.plot(V_applied_2g_poisson_x_018x1x1um, np.log10(abs(I_2g_poisson_x_018x1x1um)), 'b--', label="2 grains - Poisson extremes")

V_applied_6g_equal_doping_018x1x1um = np.delete(V_applied_6g_equal_doping_018x1x1um, 700)
I_6g_equal_doping_018x1x1um = np.delete(I_6g_equal_doping_018x1x1um, 700)
V_applied_6g_poisson_x_018x1x1um = np.delete(V_applied_6g_poisson_x_018x1x1um, 700)
I_6g_poisson_x_018x1x1um = np.delete(I_6g_poisson_x_018x1x1um, 700)

plt.plot(V_applied_6g_equal_doping_018x1x1um, np.log10(abs(I_6g_equal_doping_018x1x1um)), 'g-', label="6 grains - Equal doping")
plt.plot(V_applied_6g_poisson_x_018x1x1um, np.log10(abs(I_6g_poisson_x_018x1x1um)), 'c--', label="6 grains - Poisson extremes")
for i in range(len(V_applied_2g_equal_doping_018x1x1um)):
    print(str(f'V={V_applied_2g_equal_doping_018x1x1um[i]}, i={i}'))
# plt.title('Grain size=$0.2^{3}\mu m^{3}$')
plt.xlabel('Applied voltage (V)', fontsize=12)
plt.ylabel('log$_{10}$[Current (A)]', fontsize=12)
# plt.xlim(0, 7)
plt.ylim(-15, -5.65)
plt.legend(loc='lower right', fontsize=9)
plt.show()