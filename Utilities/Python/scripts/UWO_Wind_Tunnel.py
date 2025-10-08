#!/usr/bin/env python3
# McGrattan
# 11-21-2018
# UWO_Wind_Tunnel.py

import os
import numpy as np
import pandas as pd

outdir = '../../../out/UWO_Wind_Tunnel/'
expdir = '../../../exp/UWO_Wind_Tunnel/'

# headers
H1 = ['x', 'Cp_mean', 'Cp_rms', 'Cp_min', 'Cp_max']
H2 = ['1.0', 'NaN', 'NaN', 'NaN', 'NaN']

# read the input file (skip first header line)
input_file = os.path.join(outdir, 'UWO_inputs.csv')
C = pd.read_csv(input_file)

FDS_line_file = C.iloc[:, 0]  # first column
output_file = C.iloc[:, 1]    # second column

for k in range(len(C)):
    # read the line data (skip 2 header rows like MATLAB)
    M = pd.read_csv(os.path.join(outdir, FDS_line_file[k]), skiprows=2, header=None)

    # figure out n_sides
    n_sides = 4
    if C.iloc[k, 29] == 0: n_sides = 3
    if C.iloc[k, 20] == 0: n_sides = 2
    if C.iloc[k, 11] == 0: n_sides = 1

    # open output file for writing
    outpath = os.path.join(outdir, output_file[k])
    with open(outpath, 'w') as fid:
        fid.write(','.join(H1) + '\n')

        dist = 0.0

        for j in range(1, n_sides + 1):
            # columns in C follow the pattern in MATLAB: 3+(j-1)*9 etc.
            coord_col   = int(C.iloc[k, 2 + (j - 1) * 9]) - 1  # zero-based
            points      = int(C.iloc[k, 3 + (j - 1) * 9])
            order_index = int(C.iloc[k, 4 + (j - 1) * 9])
            x1          = float(C.iloc[k, 5 + (j - 1) * 9])
            x2          = float(C.iloc[k, 6 + (j - 1) * 9])
            mean_col    = int(C.iloc[k, 7 + (j - 1) * 9]) - 1
            rms_col     = int(C.iloc[k, 8 + (j - 1) * 9]) - 1
            min_col     = int(C.iloc[k, 9 + (j - 1) * 9]) - 1
            max_col     = int(C.iloc[k,10 + (j - 1) * 9]) - 1

            if order_index == 1:
                idx_range = range(points)
            else:
                idx_range = range(points - 1, -1, -1)

            for i in idx_range:
                x_val = (dist + M.iloc[i, coord_col] - x1) if order_index == 1 \
                        else (dist + x1 - M.iloc[i, coord_col])
                cp_mean = M.iloc[i, mean_col]
                cp_rms  = M.iloc[i, rms_col]
                cp_min  = M.iloc[i, min_col]
                cp_max  = M.iloc[i, max_col]
                fid.write(f'{x_val:6.3f},{cp_mean:6.3f},{cp_rms:6.3f},{cp_min:6.3f},{cp_max:6.3f}\n')

            # if not last side, add blank header row
            if j < n_sides:
                fid.write(','.join(H2) + '\n')

            dist += abs(x2 - x1)
