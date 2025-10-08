#!/usr/bin/env python3
# McGrattan
# 3-15-2017
# FM_Vertical_Wall_Flames.py
#
# Reads the _devc.csv file and writes output in a form appropriate for dataplot

import os
import numpy as np
import pandas as pd

outdir = '../../../out/FM_Vertical_Wall_Flames/'

nts = 32

# Read CSV skipping first two header rows
M = pd.read_csv(os.path.join(outdir, 'propylene_devc.csv'), skiprows=2, header=None)

# =============== Heat flux data ==================
H1 = ['m', 'kW/m2', 'kW/m2', 'kW/m2']
H2 = ['z', 'HF-14', 'HF-18', 'HF-20']

with open(os.path.join(outdir, 'propylene_hf.csv'), 'w') as fid:
    fid.write(','.join(H1) + '\n')
    fid.write(','.join(H2) + '\n')
    for i in range(1, 41):
        z = 0.05 * i - 0.025
        vals = M.iloc[[7, 9, 10], i].values  # MATLAB is 1-based
        fid.write(f"{z:5.3f},{vals[0]:6.2f},{vals[1]:6.2f},{vals[2]:6.2f}\n")

# =============== Temperature data ==================
H1 = ['mm', 'K', 'K', 'K', 'K', 'K']
H2 = ['x', 'T-10', 'T-12', 'T-14', 'T-18', 'T-20']

with open(os.path.join(outdir, 'propylene_tc.csv'), 'w') as fid, \
     open(os.path.join(outdir, 'propylene_tmp.csv'), 'w') as fid2:
    fid.write(','.join(H1) + '\n')
    fid.write(','.join(H2) + '\n')
    fid2.write(','.join(H1) + '\n')
    fid2.write(','.join(H2) + '\n')
    for i in range(1, 51):
        x = 3 * i - 1.5
        vals1 = M.iloc[[4, 5, 6, 8, 9], 40 + i].values + 273.15
        vals2 = M.iloc[[4, 5, 6, 8, 9], 344 + i].values + 273.15
        fid.write(f"{x:5.3f}," + ",".join([f"{v:6.2f}" for v in vals1]) + "\n")
        fid2.write(f"{x:5.3f}," + ",".join([f"{v:6.2f}" for v in vals2]) + "\n")

# =============== Soot data ==================
H1 = ['g/m2/s', 'mm', 'mm', 'mm', 'mm', 'mm']
H2 = ['mdot', '365 mm', '527 mm', '771 mm', '1022 mm', '1317 mm']

threshold = 0.0025

with open(os.path.join(outdir, 'propylene_soot.csv'), 'w') as fid:
    fid.write(','.join(H1) + '\n')
    fid.write(','.join(H2) + '\n')
    for i in range(1, nts + 1):
        delta = np.zeros(5)
        for j in range(1, 6):
            index = np.nan
            for k in range(1, 50):
                v1 = M.iloc[i, 91 + 50 * (j - 1) + k - 1]
                v2 = M.iloc[i, 91 + 50 * (j - 1) + k]
                if (v1 >= threshold) and (v2 < threshold):
                    index = k
            delta[j - 1] = index * 3 - 1.5
        mdot = 2 * i
        fid.write(f"{mdot:5.3f}," + ",".join([f"{d:5.1f}" for d in delta]) + "\n")
