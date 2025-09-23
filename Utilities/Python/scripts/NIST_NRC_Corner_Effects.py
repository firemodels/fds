#!/usr/bin/env python
"""
McGrattan
12-6-2017
NIST_NRC_Corner_Effects.py

Reads the _devc.csv file and writes output in a form appropriate for dataplot.m
"""

import pandas as pd
import numpy as np
import os

outdir = '../../../out/NIST_NRC_Corner_Effects/'

casename = [
    'corner_200_kW', 'corner_300_kW', 'corner_400_kW',
    'wall_200_kW', 'wall_300_kW', 'wall_400_kW',
    'cabinet_01', 'cabinet_02', 'cabinet_03', 'cabinet_04',
    'cabinet_05', 'cabinet_06', 'cabinet_07', 'cabinet_08',
    'cabinet_09', 'cabinet_10', 'cabinet_11', 'cabinet_12'
]

for name in casename:

    # Read CSV with pandas, skip first 2 rows to match MATLAB importdata(..., 2)
    M = pd.read_csv(os.path.join(outdir, f'{name}_devc.csv'), skiprows=2, header=None)

    # header info
    H1 = ['s', 'C', 'C', 'C']
    H2 = ['Time', 'Lower', 'Middle', 'Upper']

    with open(os.path.join(outdir, f'{name}_Plume.csv'), 'w', newline='') as fid:
        fid.write(','.join(H1) + '\n')
        fid.write(','.join(H2) + '\n')

        n_times = len(M)
        for i in range(n_times):
            # smoothing window ±5 rows
            start_idx = max(0, i - 5)
            end_idx   = min(n_times, i + 6)  # Python slicing excludes last index

            # MATLAB column indices: 60:88 → Python 59:88, 31:59 → 30:59, 2:30 → 1:30
            low_smooth  = M.iloc[start_idx:end_idx, 59:88].mean(axis=0)
            mid_smooth  = M.iloc[start_idx:end_idx, 30:59].mean(axis=0)
            high_smooth = M.iloc[start_idx:end_idx, 1:30].mean(axis=0)

            low  = low_smooth.max()
            mid  = mid_smooth.max()
            high = high_smooth.max()

            fid.write(f'{int(M.iloc[i,0]):4d},{low:5.1f},{mid:5.1f},{high:5.1f}\n')

    del M
