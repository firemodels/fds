
import os
import csv
import numpy as np
import pandas as pd

# McGrattan
# 07-02-2025
# UWO_Wind Tunnel.py
#
# Reads and transposes line data for UWO Wind Tunnel cases

outdir = '../../../out/UWO_Wind_Tunnel/'
expdir = '../../../exp/UWO_Wind_Tunnel/'

H = [['x', 'Cp_mean', 'Cp_rms', 'Cp_min', 'Cp_max'],
     ['1.0', 'NaN', 'NaN', 'NaN', 'NaN']]

# Read input CSV file
input_file_path = os.path.join(outdir, 'UWO_inputs.csv')
with open(input_file_path, 'r') as input_file:
    csv_reader = csv.reader(input_file)
    next(csv_reader)  # Skip header line
    rows = list(csv_reader)

# Parse the CSV data into columns (equivalent to textscan)
C = [[] for _ in range(38)]  # 38 columns based on the textscan format string

for row in rows:
    C[0].append(row[0])   # %s - FDS_line_file
    C[1].append(row[1])   # %s - output_file
    for i in range(2, 38):  # remaining columns are integers or floats
        if i in [5, 6, 14, 15, 23, 24, 32, 33]:  # %f columns
            C[i].append(float(row[i]))
        else:  # %d8 columns
            C[i].append(int(row[i]))

FDS_line_file = C[0]
output_file = C[1]

for k in range(len(C[0])):
    
    # Import data (equivalent to importdata)
    data_file_path = os.path.join(outdir, C[0][k])
    M_data = pd.read_csv(data_file_path, skiprows=2, header=None).values
    
    n_sides = 4
    if C[29][k] == 0:  # C{30}(k) in MATLAB (0-indexed in Python)
        n_sides = 3
    if C[20][k] == 0:  # C{21}(k) in MATLAB
        n_sides = 2
    if C[11][k] == 0:  # C{12}(k) in MATLAB
        n_sides = 1
    
    # Open output file for writing
    output_file_path = os.path.join(outdir, output_file[k])
    with open(output_file_path, 'w', newline='') as fid:
        # Write header
        fid.write(','.join(H[0]) + '\n')
        
        dist = 0.0
        
        for j in range(n_sides):  # j from 0 to n_sides-1 in Python
            
            coord_col = C[2 + j*9][k]      # C{ 3+(j-1)*9}(k) in MATLAB
            points = C[3 + j*9][k]         # C{ 4+(j-1)*9}(k) in MATLAB
            order_index = C[4 + j*9][k]    # C{ 5+(j-1)*9}(k) in MATLAB
            x1 = C[5 + j*9][k]             # C{ 6+(j-1)*9}(k) in MATLAB
            x2 = C[6 + j*9][k]             # C{ 7+(j-1)*9}(k) in MATLAB
            mean_col = C[7 + j*9][k]       # C{ 8+(j-1)*9}(k) in MATLAB
            rms_col = C[8 + j*9][k]        # C{ 9+(j-1)*9}(k) in MATLAB
            min_col = C[9 + j*9][k]        # C{10+(j-1)*9}(k) in MATLAB
            max_col = C[10 + j*9][k]       # C{11+(j-1)*9}(k) in MATLAB
            
            # Adjust column indices for 0-based indexing
            coord_col -= 1
            mean_col -= 1
            rms_col -= 1
            min_col -= 1
            max_col -= 1
            
            if order_index == 1:
                for i in range(points):  # i from 0 to points-1 in Python
                    line_data = [
                        dist + M_data[i, coord_col] - x1,
                        M_data[i, mean_col],
                        M_data[i, rms_col],
                        M_data[i, min_col],
                        M_data[i, max_col]
                    ]
                    fid.write('{:6.3f},{:6.3f},{:6.3f},{:6.3f},{:6.3f}\n'.format(*line_data))
            else:
                for i in range(points-1, -1, -1):  # i from points-1 down to 0
                    line_data = [
                        dist + x1 - M_data[i, coord_col],
                        M_data[i, mean_col],
                        M_data[i, rms_col],
                        M_data[i, min_col],
                        M_data[i, max_col]
                    ]
                    fid.write('{:6.3f},{:6.3f},{:6.3f},{:6.3f},{:6.3f}\n'.format(*line_data))
            
            if j < n_sides - 1:  # j+1 < n_sides in MATLAB equivalent
                fid.write(','.join(H[1]) + '\n')
            
            dist = dist + abs(x2 - x1)
    
    # Clear M equivalent (Python garbage collection handles this automatically)
    del M_data


