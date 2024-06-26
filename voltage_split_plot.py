import matplotlib.pyplot as plt
import numpy as np
import os
import csv

path_here = os.path.dirname(os.path.realpath(__file__))

doping = [1e16, 1.2e16, 1.4e16, 1.6e16, 1.8e16, 2e16]

VA = [-0.9,-0.6,-0.3,0.3,0.6,0.9]

doping_file = path_here + "\\voltage_split_data_doping.csv"

alpha = np.zeros((len(doping), len(VA)))

counter = 0

stop_outer_loop = False

with open(doping_file, 'r') as csvfile:
    csvreader = csv.reader(csvfile, delimiter=';')
    next(csvreader)
    for i, row in enumerate(csvreader):
        
        current_series = []
        # Iterate over the next 6 rows to extract values and append them to the current series
        for j in range(6):
            # Append the value of the ninth column to the current series
            try:
                current_series.append(eval(row[8]))
                # Move to the next row
                if j != 5:
                    row = next(csvreader)
            except:
                stop_outer_loop = True
                break
        if stop_outer_loop:
            break
        # Append the current series to the series_values array
        current_series =np.array(current_series).T
        print(current_series.shape)
        alpha[:,counter] = current_series
        counter+=1

# Print the array containing the values from each series of 6 rows
print(alpha)


for v in range(len(VA)):
    print(f'VA={VA[v]}, alpha = {alpha[:,v]}')
    plt.plot(doping, alpha[:,v], label='$V_{A}$'+f'={VA[v]}V')
plt.xlabel('Doping in RHS', fontsize=16)
plt.ylabel('$\\alpha$', fontsize=16)
# plt.legend(fontsize=16)
plt.figlegend(loc='upper left', bbox_to_anchor=(0.08, 1), ncol=3, fontsize=15)
# plt.figlegend(loc="center right", bbox_to_anchor=(1.05, 0.5), fontsize=16)
plt.tight_layout()


# consider different doping
doping2 = np.concatenate((np.arange(1e16, 2.2e16, 2e15), np.arange(1.2e16, 2.2e16, 2e15)))
VA2 = np.concatenate((np.arange(-0.9,0,0.1), np.arange(0.1,1,0.1)))

voltage_file =  path_here + "\\voltage_split_data_edited.csv"

alpha2 = np.zeros((len(VA2), len(doping2)))

counter = 0

stop_outer_loop = False

with open(voltage_file, 'r') as csvfile:
    csvreader = csv.reader(csvfile, delimiter=';')
    next(csvreader)
    next(csvreader)
    for i, row in enumerate(csvreader):
        
        current_series = []
        # Iterate over the next 6 rows to extract values and append them to the current series
        for j in range(18):
            # Append the value of the ninth column to the current series
            current_series.append(eval(row[8]))
            # Move to the next row
            try:
                row = next(csvreader)
            except:
                print("Final row")
            if counter == 12:
                stop_outer_loop = True
                break
        if stop_outer_loop:
            break
        # Append the current series to the series_values array
        current_series =np.array(current_series).T
        print(current_series.shape)
        alpha2[:,counter] = current_series
        print(f'ND={doping2[counter]}')
        counter+=1

names = ["ND_{1}=$1\\times 10^{16}$, ND_{2}=$1\\times 10^{16}$",
         "ND_{1}=$1\\times 10^{16}$, ND_{2}=$1.6\\times 10^{16}$",
         "ND_{1}=$1\\times 10^{16}$, ND_{2}=$2\\times 10^{16}$",
         "ND_{1}=$1.6\\times 10^{16}$, ND_{2}=$1\\times 10^{16}$",
         "ND_{1}=$2\\times 10^{16}$, ND_{2}=$1\\times 10^{16}$"]

names = ["$1\\times 10^{16}$ cm$^{-3}$","$1.6\\times 10^{16}$ cm$^{-3}$","$2\\times 10^{16}$ cm$^{-3}$"]

plt.figure()
counter=0
for ND in [0,3,5]:
    print(f'ND2={doping2[ND]}, alpha = {alpha2[:,ND]}')
    plt.plot(VA2, alpha2[:,ND], label='N$_{D}$'+f'={names[counter]}')
    counter+=1
plt.xlabel('Applied voltage', fontsize=16)
plt.ylabel('$\\alpha$', fontsize=16)
plt.legend(loc='upper right', fontsize=12)
plt.show()