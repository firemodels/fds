#!/usr/bin/env python3
# McGrattan
# 2-26-2018
# umd_line_burner_process.py
#
# Read and process FDS output files for UMD Line Burner

import os
import numpy as np
import pandas as pd

Lf_dt = 10

outdir = '../../../out/UMD_Line_Burner/'

fuel_name = ['methane', 'propane']
res_name = ['1p25cm', 'p625cm', 'p3125cm']

for i_fuel in range(2):

    for fds_resolution in range(3):

        dev_file = os.path.join(outdir, f"{fuel_name[i_fuel]}_dx_{res_name[fds_resolution]}_devc.csv")
        hrr_file = os.path.join(outdir, f"{fuel_name[i_fuel]}_dx_{res_name[fds_resolution]}_hrr.csv")

        # Read CSVs, skipping first two header lines
        DEV = pd.read_csv(dev_file, skiprows=2, header=None)
        HRR = pd.read_csv(hrr_file, skiprows=2, header=None)

        # Read headers separately (the first row after skiprows=1)
        with open(dev_file, 'r') as f:
            header_lines = [next(f) for _ in range(2)]
        dev_headers = header_lines[1].strip().split(',')

        with open(hrr_file, 'r') as f:
            header_lines = [next(f) for _ in range(2)]
        hrr_headers = header_lines[1].strip().split(',')

        # Locate columns (match MATLAB’s strcmp behavior exactly)
        Time_idx = dev_headers.index('Time')
        XO2_idx = dev_headers.index('XO2')
        qrad1_idx = dev_headers.index('qrad1')
        qrad2_idx = dev_headers.index('qrad2')
        Lf_idx = dev_headers.index('L_F')

        HRR_idx = hrr_headers.index('HRR')
        Qrad_idx = hrr_headers.index('Q_RADI')

        Time_FDS = DEV.iloc[:, Time_idx].values
        XO2_FDS = DEV.iloc[:, XO2_idx].values
        Qdot_FDS = HRR.iloc[:, HRR_idx].values
        Qrad_FDS = HRR.iloc[:, Qrad_idx].values
        q_R_FDS = 0.5 * (DEV.iloc[:, qrad1_idx].values + DEV.iloc[:, qrad2_idx].values)
        Lf_FDS = DEV.iloc[:, Lf_idx].values.copy()

        ntp = len(Time_FDS)

        # Moving-average Lf_FDS with Lf_dt window
        Lf_tmp = Lf_FDS.copy()
        for n in range(ntp):
            mask = Time_FDS > (Time_FDS[n] - Lf_dt)
            # indices from first True to n inclusive
            idxs = np.where(mask)[0]
            idxs = idxs[idxs <= n]
            Lf_FDS[n] = np.mean(Lf_tmp[idxs])

        # Write output file
        out_file = os.path.join(outdir, f"{fuel_name[i_fuel]}_dx_{res_name[fds_resolution]}.csv")
        with open(out_file, 'w') as fid:
            fid.write('XO2,eta,Chi_R,Lf,q_R\n')
            for ii in range(ntp):
                eta = Qdot_FDS[ii] / 50.0
                Chi_R = max(0, -Qrad_FDS[ii] / max(0.001, Qdot_FDS[ii]))
                fid.write(f"{XO2_FDS[ii]:5.3f},{eta:6.2f},{Chi_R:6.2f},{Lf_FDS[ii]:6.2f},{q_R_FDS[ii]:6.2f}\n")
