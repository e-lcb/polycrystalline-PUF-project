import matplotlib.pyplot as plt
import numpy as np
import ltspice
import os

path_here = os.path.dirname(os.path.realpath(__file__))

#MAKE SURE ALL RAW FILES ARE RUN IN THE SAME RANGE

fp_15x15_point1="C:\\Users\\elc20\\OneDrive\\Documents\\University\\Year 4\\Project\\2D Model\\Report things 2 Random length and doping\\Si 15x15 ND 1e16 L 15um T 300K NT 4e11 Point Contacts.raw"
fp_15x15_point2="C:\\Users\\elc20\\OneDrive\\Documents\\University\\Year 4\\Project\\2D Model\\Report things 2 Random length and doping\\Si 15x15 ND 1e16 L 15um T 300K NT 4e11 Point Contacts 2.raw"
fp_15x15_point3="C:\\Users\\elc20\\OneDrive\\Documents\\University\\Year 4\\Project\\2D Model\\Report things 2 Random length and doping\\Si 15x15 ND 1e16 L 15um T 300K NT 4e11 Point Contacts 3.raw"
fp_15x15_point4="C:\\Users\\elc20\\OneDrive\\Documents\\University\\Year 4\\Project\\2D Model\\Report things 2 Random length and doping\\Si 15x15 ND 1e16 L 15um T 300K NT 4e11 Point Contacts 4.raw"
fp_15x15_point5="C:\\Users\\elc20\\OneDrive\\Documents\\University\\Year 4\\Project\\2D Model\\Report things 2 Random length and doping\\Si 15x15 ND 1e16 L 15um T 300K NT 4e11 Point Contacts 5.raw"
fp_15x15_point6="C:\\Users\\elc20\\OneDrive\\Documents\\University\\Year 4\\Project\\2D Model\\Report things 2 Random length and doping\\Si 15x15 ND 1e16 L 15um T 300K NT 4e11 Point Contacts 6.raw"
fp_15x15_point7="C:\\Users\\elc20\\OneDrive\\Documents\\University\\Year 4\\Project\\2D Model\\Report things 2 Random length and doping\\Si 15x15 ND 1e16 L 15um T 300K NT 4e11 Point Contacts 7.raw"
fp_15x15_point8="C:\\Users\\elc20\\OneDrive\\Documents\\University\\Year 4\\Project\\2D Model\\Report things 2 Random length and doping\\Si 15x15 ND 1e16 L 15um T 300K NT 4e11 Point Contacts 8.raw"
fp_15x15_point9="C:\\Users\\elc20\\OneDrive\\Documents\\University\\Year 4\\Project\\2D Model\\Report things 2 Random length and doping\\Si 15x15 ND 1e16 L 15um T 300K NT 4e11 Point Contacts 9.raw"
fp_15x15_point10="C:\\Users\\elc20\\OneDrive\\Documents\\University\\Year 4\\Project\\2D Model\\Report things 2 Random length and doping\\Si 15x15 ND 1e16 L 15um T 300K NT 4e11 Point Contacts 10.raw"
fp_15x15_point11="C:\\Users\\elc20\\OneDrive\\Documents\\University\\Year 4\\Project\\2D Model\\Report things 2 Random length and doping\\Si 15x15 ND 1e16 L 15um T 300K NT 4e11 Point Contacts 11.raw"
fp_15x15_point12="C:\\Users\\elc20\\OneDrive\\Documents\\University\\Year 4\\Project\\2D Model\\Report things 2 Random length and doping\\Si 15x15 ND 1e16 L 15um T 300K NT 4e11 Point Contacts 12.raw"
fp_15x15_point13="C:\\Users\\elc20\\OneDrive\\Documents\\University\\Year 4\\Project\\2D Model\\Report things 2 Random length and doping\\Si 15x15 ND 1e16 L 15um T 300K NT 4e11 Point Contacts 13.raw"
fp_15x15_point14="C:\\Users\\elc20\\OneDrive\\Documents\\University\\Year 4\\Project\\2D Model\\Report things 2 Random length and doping\\Si 15x15 ND 1e16 L 15um T 300K NT 4e11 Point Contacts 14.raw"

fp_50x50_point="C:\\Users\\elc20\\OneDrive\\Documents\\University\\Year 4\\Project\\2D Model\\Report things 2 Random length and doping\\Si 50x50 ND 1e16 L 15um T 300K NT 4e11 Point Contacts.raw"


fp_15x15_1_in_0__7="C:\\Users\\elc20\\OneDrive\\Documents\\University\\Year 4\\Project\\2D Model\\Report things 2 Random length and doping\\Si 15x15 ND 1e16 L 15um T 300K NT 4e11 [0,7] Contacts 1.raw"
fp_15x15_1_in_14="C:\\Users\\elc20\\OneDrive\\Documents\\University\\Year 4\\Project\\2D Model\\Report things 2 Random length and doping\\Si 15x15 ND 1e16 L 15um T 300K NT 4e11 in 14 Contacts 1.raw"
fp_15x15_1_in_7="C:\\Users\\elc20\\OneDrive\\Documents\\University\\Year 4\\Project\\2D Model\\Report things 2 Random length and doping\\Si 15x15 ND 1e16 L 15um T 300K NT 4e11 in 7 Contacts 1.raw"
fp_15x15_1_in_3__8="C:\\Users\\elc20\\OneDrive\\Documents\\University\\Year 4\\Project\\2D Model\\Report things 2 Random length and doping\\Si 15x15 ND 1e16 L 15um T 300K NT 4e11 in [3,8] Contacts 1.raw"
fp_15x15_1_in_0__14="C:\\Users\\elc20\\OneDrive\\Documents\\University\\Year 4\\Project\\2D Model\\Report things 2 Random length and doping\\Si 15x15 ND 1e16 L 15um T 300K NT 4e11 in [0,14] Contacts 1.raw"
fp_15x15_1_swapped_point1="C:\\Users\\elc20\\OneDrive\\Documents\\University\\Year 4\\Project\\2D Model\\Report things 2 Random length and doping\\Si 15x15 ND 1e16 L 15um T 300K NT 4e11 swapped Point Contacts 1.raw"

fp_15x15_1_in_electrode="C:\\Users\\elc20\\OneDrive\\Documents\\University\\Year 4\\Project\\2D Model\\Report things 2 Random length and doping\\Si 15x15 ND 1e16 L 15um T 300K NT 4e11 in electrode Contacts 1.raw"

fp_15x15_in_0_5__10_15_out_0__14="C:\\Users\\elc20\\OneDrive\\Documents\\University\\Year 4\\Project\\2D Model\\Report things 2 Random length and doping\\Si 15x15 ND 1e16 L 15um T 300K NT 4e11 Point Contacts MIMG (0-5,10-15) (0,14).raw"

fp_15x15_point1_GaAs="C:\\Users\\elc20\\OneDrive\\Documents\\University\\Year 4\\Project\\2D Model\\Report things 3 Different Material\\GaAs 15x15 ND 1e16 L 15um T 300K NT 4e11 Point Contacts.raw"

fp_15x15_point1_358K = "C:\\Users\\elc20\\OneDrive\\Documents\\University\\Year 4\\Project\\2D Model\\Report things 4 Different temperatures\\Si 15x15 ND 1e16 L 15um T 358K NT 4e11 Point Contacts.raw"
fp_15x15_point1_325K = "C:\\Users\\elc20\\OneDrive\\Documents\\University\\Year 4\\Project\\2D Model\\Report things 4 Different temperatures\\Si 15x15 ND 1e16 L 15um T 325K NT 4e11 Point Contacts.raw"
fp_15x15_point1_343K = "C:\\Users\\elc20\\OneDrive\\Documents\\University\\Year 4\\Project\\2D Model\\Report things 4 Different temperatures\\Si 15x15 ND 1e16 L 15um T 343K NT 4e11 Point Contacts.raw"
fp_15x15_point1_273K = "C:\\Users\\elc20\\OneDrive\\Documents\\University\\Year 4\\Project\\2D Model\\Report things 4 Different temperatures\\Si 15x15 ND 1e16 L 15um T 273K NT 4e11 Point Contacts.raw"

fp_15x15_in_0_4__5_9__10_15_out_14="C:\\Users\\elc20\\OneDrive\\Documents\\University\\Year 4\\Project\\2D Model\\Report things 2 Random length and doping\\Si 15x15 ND 1e16 L 15um T 300K NT 4e11 Point Contacts MIMG (0-4,5-9,10-15) (14).raw"
fp_15x15_3_in_0_4__5_9__10_15_out_14="C:\\Users\\elc20\\OneDrive\\Documents\\University\\Year 4\\Project\\2D Model\\Report things 2 Random length and doping\\Si 15x15 ND 1e16 L 15um T 300K NT 4e11 MIMG (0-4,5-9,10-15) (14) 3.raw"

# fp_15x15_3_in_5bits_out_14="C:\\Users\\elc20\\OneDrive\\Documents\\University\\Year 4\\Project\\2D Model\\Report things 2 Random length and doping\\Si 15x15 ND 1e16 L 15um T 273K NT 4e11 MIMG 5 in (14) 3.raw"

# fp_15x15_1_in_5bits_out_14="C:\\Users\\elc20\\OneDrive\\Documents\\University\\Year 4\\Project\\2D Model\\Report things 2 Random length and doping\\Si 15x15 ND 1e16 L 15um T 300K NT 4e11 MIMG 5 in (14) 1.raw"
fp_15x15_2_in_5bits_out_14="C:\\Users\\elc20\\OneDrive\\Documents\\University\\Year 4\\Project\\2D Model\\Report things 2 Random length and doping\\Si 15x15 ND 1e16 L 15um T 300K NT 4e11 MIMG 5 in (14) 2.raw"
fp_15x15_5_in_5bits_out_14="C:\\Users\\elc20\\OneDrive\\Documents\\University\\Year 4\\Project\\2D Model\\Report things 2 Random length and doping\\Si 15x15 ND 1e16 L 15um T 300K NT 4e11 MIMG 5 in (14) 5.raw"


# raw_files = [fp_15x15_point1,fp_15x15_point2,fp_15x15_in_0_5__10_15_out_0__14,fp_15x15_in_0_4__5_9__10_15_out_14, fp_15x15_point1_GaAs]
raw_files = [fp_15x15_point1,fp_15x15_point2,fp_15x15_point3, fp_15x15_point4, fp_15x15_point5, fp_15x15_point6, fp_15x15_point7, fp_15x15_point8, fp_15x15_point9, fp_15x15_point10, fp_15x15_point11, fp_15x15_point12, fp_15x15_point13, fp_15x15_point14] # random points
# raw_files = [fp_15x15_point1,fp_15x15_point2,fp_15x15_point3, fp_15x15_point4, fp_15x15_point5] # random points
# raw_files = [fp_15x15_point1_273K, fp_15x15_point1, fp_15x15_point1_325K, fp_15x15_point1_343K, fp_15x15_point1_358K] #different temperatures
# raw_files = [fp_15x15_in_0_4__5_9__10_15_out_14, fp_15x15_3_in_0_4__5_9__10_15_out_14]

# raw_files = [fp_15x15_point1, fp_15x15_1_in_14, fp_15x15_1_in_7, fp_15x15_1_in_3__8, fp_15x15_1_in_0__14, fp_15x15_1_swapped_point1] # different inputs

# raw_files = [fp_15x15_point1, fp_15x15_1_in_7, fp_15x15_1_in_14]

# raw_files = [fp_15x15_point1,fp_15x15_1_in_electrode]

# raw_files = [fp_15x15_2_in_5bits_out_14]
# raw_files=[fp_15x15_point1]

raw_files = [fp_50x50_point]

m = [15,15,15,15,15]
n = [15,15,15,15,15]

n=[50]
m=[50]
# VA = np.linspace(0.9,1.1,21)
VA=[15,2,5,10,15]
graph_V = 0 #V index of VA you want to see
filenames=list(range(len(raw_files)))
filename_graph=list(range(len(raw_files)))
# filename_graph=["bottom", "middle", "top"]

save_dir = "fp_15x15_2_in_5bits_out_14"
save_file = None #accidentally overwrote 00001 fp_15x15_2_in_5bits_out_14



def nearest_voltage(raw_file, VA):
    filepath = raw_file
    filesize = os.stat(filepath).st_size
    print(str(f'filesize = {filesize}'))
    LT = ltspice.Ltspice(filepath)
    LT.parse()
    V_source = LT.get_data('V1')
    
    indx = np.searchsorted(V_source, VA) #the voltage may not exactly match the desired applied voltage so have to find the closest one
    
    if abs(V_source[indx]-VA) < abs(V_source[indx-1]-VA):
        indx = indx
    else:
        indx = indx - 1
    VA_actual = V_source[indx] #this is the closest value from the simulation
    # np.savetxt("C:\\Users\\elc20\\OneDrive\\Documents\\University\\Year 4\\Project\\2D Model\\V_source_test.txt", V_source)
    print(str(f'index = {indx}'))
    return indx, VA_actual

def read_file(raw_file):
    filepath = raw_file
    name = ltspice.Ltspice(filepath)
    name.parse()
    return name

indx = []

for i in range(len(VA)):
    indx.append(nearest_voltage(VA=VA[i],raw_file=raw_files[0])[0])

VA_actual=[]

for j in range(len(raw_files)):
    filenames[j]=read_file(raw_files[j])
    for i in range(len(indx)):
        VA_actual.append(filenames[j].get_data('V1')[indx[i]])

print(f'VA_actual for all sims: {VA_actual}')

edge = 3*m[0]-3
V_edge = np.zeros((n[0]-1, len(VA),len(raw_files)))
I_edge = np.zeros((n[0]-1, len(VA),len(raw_files)))
f=0

for name in filenames:
    for i in range(n[0]-1):
        for ind in range(len(indx)):
            # print(str(f'(Rx{edge}y{i}), index={indx[ind]}'))
            try:
                V_edge[i,ind,f] = name.get_data(str(f'V(x{edge:02.0f}y{i:02.0f})'))[indx[ind]]
            except:
                try:
                    V_edge[i,ind,f] = name.get_data(str(f'V(NV5)'))[indx[ind]]
                except:
                    V_edge[i,ind,f] = None
                
            try:
                I_edge[i,ind,f] = name.get_data(str(f'I(Rx{edge-1}y{i})'))[indx[ind]]
            except:
                I_edge[i,ind,f] = None
    f+=1

f = 0
counter = 0
bottom = n[0] -1
V_bottom = np.zeros((len(range(0,edge+1,3)), len(VA),len(raw_files)))
I_bottom = np.zeros((len(range(0,edge+1,3)), len(VA),len(raw_files)))
for name in filenames:
    counter=0
    for i in range(0,edge+1,3):
        for ind in range(len(indx)):
            # print(str(f'(Rx{i}y{bottom}), index={indx[ind]}, counter={counter}'))
            try:
                V_bottom[counter,ind,f] = name.get_data(str(f'V(x{i:02.0f}y{bottom:02.0f})'))[indx[ind]]
            except:
                V_bottom[counter,ind,f] = None

            try:
                I_bottom[counter,ind,f] = name.get_data(str(f'I(Rx{i}y{bottom})'))[indx[ind]]
            except:
                I_bottom[counter,ind,f] = None
        counter+=1
    f+=1

bottom_length=[]
for i in range(len(range(0,edge+1,3))):
    bottom_length.append(i)

if save_file != None:
    path = path_here + "\\Report things 2 Random length and doping\\"+ save_dir
    if not os.path.exists(path):
        os.mkdir(path)
    path += "\\"+save_file+".csv"
    np.savetxt(path, V_edge[:,graph_V,0], delimiter=',')
    

plt.figure(1)
#Different sims (constant voltage)
for f in range(len(raw_files)):
    # print(str(V_edge[:,graph_V,f]))
    plt.plot(range(n[0]-1), V_edge[:,graph_V,f], label=f+1)
plt.xlabel("Grain number along right edge.",fontsize=16)
plt.ylabel("Voltage (V)", fontsize=16)
plt.legend()

plt.figure(2)
for f in range(len(raw_files)):
    # print(str(I_edge[:,graph_V,f]))
    plt.plot(range(n[0]-1), I_edge[:,graph_V,f], label=f+1)
plt.xlabel("Grain number along right edge.",fontsize=16)
plt.ylabel("Current (A)", fontsize=16)
plt.legend()

plt.figure(3)
for f in range(len(raw_files)):
    plt.plot(bottom_length, V_bottom[:,graph_V,f], label=f+1)
plt.xlabel("Grain number along top edge.",fontsize=16)
plt.ylabel("Voltage (V)", fontsize=16)
plt.legend()

plt.figure(4)
for f in range(len(raw_files)):
    plt.plot(bottom_length, I_bottom[:,graph_V,f], label=f+1)
plt.xlabel("Grain number along top edge.",fontsize=16)
plt.ylabel("Current (A)", fontsize=16)
plt.legend()

#Different voltages (same sim)
plt.figure(5)
for v in range(len(VA)):
    plt.plot(range(n[0]-1), V_edge[:,v,0],label=str(f'VA={VA[v]}V'))
plt.xlabel("Grain number along right edge.",fontsize=16)
plt.ylabel("Voltage (V)", fontsize=16)
plt.legend(fontsize=16)

plt.figure(6)
for v in range(len(VA)):
    plt.plot(bottom_length, V_bottom[:,v,0],label=str(f'VA={VA[v]}V'))
plt.xlabel("Grain number along bottom edge.",fontsize=16)
plt.ylabel("Voltage (V)", fontsize=16)
plt.legend(fontsize=16)

plt.show()


error = np.zeros((n[0]-1, len(VA),len(raw_files)))
error_current = np.zeros((n[0]-1, len(VA),len(raw_files)))

# for f in range(len(filenames)):
#     for i in range(n[0]-1):
#         for v in range(len(VA)):
#             error[i,v,f] = abs(V_edge[i,v,f]-V_edge[i,10,f])*100/abs(V_edge[i,10,f])

# plt.figure(7)

names = ["top", "middle", "bottom"]
# counter = 0
# for i in [0,7,13]:
#     plt.plot(VA, error[i,:,0], label=names[counter])
#     counter+=1
# plt.xlabel("Applied voltage (V)",fontsize=16)
# plt.ylabel("% error", fontsize=16)
# plt.xticks([0.9,0.95,1,1.05,1.1])
# plt.legend(fontsize=16)

# #TEMPERATURE PLOTS
# temperature = [273,300,325,343,358]

# for f in range(len(raw_files)):
#     for i in range(n[0]-1):
#         error[i,0,f] = abs(V_edge[i,0,f]-V_edge[i,0,1])*100/abs(V_edge[i,0,1])
#         error_current[i,0,f] = abs(I_edge[i,0,f]-I_edge[i,0,1])*100/abs(I_edge[i,0,1])


# plt.figure(8)
# counter=0
# for i in [0,7,13]:
#     plt.plot(temperature, error[i,0,:], label=names[counter])
#     plt.plot(temperature, error[i,0,:], ".k")
#     counter+=1
# plt.xlabel("Temperature",fontsize=16)
# plt.ylabel("% error in voltage", fontsize=16)
# plt.xticks(temperature)
# plt.legend(fontsize=16)

# plt.figure(9)
# counter=0
# for i in [0,7,13]:
#     plt.plot(temperature, error_current[i,0,:], label=names[counter])
#     plt.plot(temperature, error_current[i,0,:], ".k")
#     counter+=1
# plt.xlabel("Temperature",fontsize=16)
# plt.ylabel("% error in current", fontsize=16)
# plt.xticks(temperature)
# plt.legend(fontsize=16)


#Use quantisation on the edge for a couple of middle grains
# Define parameters
min_voltage = 0  # Minimum voltage
max_voltage = 0.2  # Maximum voltage
resolution = 0.05  # Voltage resolution

# Calculate number of quantization levels
num_levels = int((max_voltage - min_voltage) / resolution) + 1

# Calculate number of bits needed to represent all levels
num_bits = len(bin(num_levels - 1)) - 2  # Subtract 2 for '0b' prefix

# Generate binary codes for each level
binary_codes = {round(min_voltage + i * resolution, 2): bin(i)[2:].zfill(num_bits) for i in range(num_levels)}

# Print the binary codes for each voltage level
# for voltage, code in binary_codes.items():
#     print(f"Voltage: {voltage}V, Binary Code: {code}")


#Voltage to quantise
# edge_points=[5, 7, 9]
edge_points = list(range(14))
bottom_points = bottom_length[:-1]

V_quant = V_edge[edge_points, 0, :]
V_quant_bottom = V_bottom[bottom_points, 0, :]

print(f'V_quant={V_quant}')
print(f'V_quant_bottom={V_quant_bottom}')



binary_code=[]

binary_array = np.zeros((len(edge_points), len(raw_files)), dtype='S' + str(num_bits))
binary_array_bottom = np.zeros((len(bottom_points), len(raw_files)), dtype='S' + str(num_bits))


# Calculate the quantization level

for f in range(len(raw_files)):
    for pos in range(len(edge_points)):
        quantization_level = int((V_quant[pos,f] - min_voltage) / resolution + 0.5)

        # Retrieve the binary code from the lookup table
        binary_code.append(binary_codes.get(round(min_voltage + quantization_level * resolution, 2), "Not found"))
        
        print(f"The binary representation of {V_quant[pos,f]}V is: {binary_code[pos]}")

    binary_array[:,f]=binary_code
    binary_code=[]

for f in range(len(raw_files)):
    for pos in range(len(bottom_points)):
        quantization_level = int((V_quant_bottom[pos,f] - min_voltage) / resolution + 0.5)

        # Retrieve the binary code from the lookup table
        binary_code.append(binary_codes.get(round(min_voltage + quantization_level * resolution, 2), "Not found"))
        
        print(f"The binary representation of {V_quant_bottom[pos,f]}V is: {binary_code[pos]} (bottom)")

    binary_array_bottom[:,f]=binary_code
    binary_code=[]

print(binary_array)

joined_array=np.zeros(len(raw_files), dtype='S' + str(num_bits*len(raw_files)))
joined_array_bottom=np.zeros(len(raw_files), dtype='S' + str(num_bits*len(raw_files)))
joined_array_both=np.zeros(len(raw_files), dtype='S' + str(num_bits*len(raw_files)))


for i in range(np.size(binary_array,1)):
    joined_array[i] = b''.join(binary_array[:,i])

for i in range(np.size(binary_array_bottom,1)):
    joined_array_bottom[i] = b''.join(binary_array_bottom[:,i])

print(f'joined_array={joined_array}')
print(f'joined_array_bottom={joined_array_bottom}')

# for i in range(len(raw_files)):
#     joined_array[i] = b''.join(joined_array[i],joined_array_bottom[i])

joined_array_both = [b1 + b2 for b1, b2 in zip(joined_array, joined_array_bottom)]

print(f'joined_array_both={joined_array_both}')

# Function to calculate Hamming distance between two byte strings
def hamming_distance(byte_str1, byte_str2):
    # Ensure both byte strings have the same length
    if len(byte_str1) != len(byte_str2):
        raise ValueError("Byte strings must have the same length")

    # Calculate the Hamming distance
    distance = sum(c1 != c2 for c1, c2 in zip(byte_str1, byte_str2))
    return distance


#REMOVE THIS TO GO BACK TO JUST RIGHT EDGE
joined_array=joined_array_both

# Calculate Hamming distance between each pair of byte strings
num_strings = len(joined_array)
hamming_distances = np.zeros((num_strings, num_strings), dtype=int)

for i in range(num_strings):
    for j in range(num_strings):
        hamming_distances[i, j] = hamming_distance(joined_array[i], joined_array[j])

# Print the Hamming distances
print("Hamming distances between byte strings:")
print(hamming_distances)

num_bits_total = len(joined_array[0])
print(f'total num bits = {num_bits_total}')

k=len(raw_files)
row_sums = np.sum(hamming_distances, axis=1)
total_sum = np.sum(row_sums)/(num_bits_total)*100
uniqueness = 2/(k*(k-1))*total_sum

reliability = 100 - row_sums/(num_bits_total)*100/(k-1) #k-1 is the number of other options being compared

print(f'Uniqueness={uniqueness}%, k={k}, m={num_bits*len(edge_points)}')
print(f'Reliability={reliability}%')

print(f'Mean reliability: {np.mean(reliability)}%')
print(f'Intra HD: {100-np.mean(reliability)}%')

def calculate_uniformity(byte_strings):
    uniformities = []
    for byte_string in byte_strings:
        # Convert byte string to binary string representation
        # binary_string = ''.join(format(byte, '08b') for byte in byte_string)
        binary_string=byte_string.decode('utf-8')
        # print(f'byte string = {binary_string}')
        # Calculate uniformity using the binary string
        num_zeros = binary_string.count('0')
        num_ones = binary_string.count('1')
        total_bits = len(binary_string)
        # print(f'num zeros={num_zeros}, num_ones={num_ones}')
        uniformity = min(num_zeros, num_ones) / total_bits *100
        uniformities.append(uniformity)
    return uniformities

# Calculate uniformity for each byte string
uniformities = calculate_uniformity(joined_array)

# Print the uniformities
print("Uniformity of each byte string:")
for i, uniformity in enumerate(uniformities, 1):
    print(f"Byte string {i}: {uniformity}")

average_uniformity = np.mean(uniformities)

print(f'Mean uniformity: {average_uniformity}%')

# plt.show()