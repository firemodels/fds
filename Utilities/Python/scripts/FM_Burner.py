"""
FM_Burner.py
Python translation of FM_Burner.m (McGrattan, 7-30-2018)

Reads and processes FDS output files for FM_Burner cases.
"""

import os
import numpy as np
import pandas as pd

outdir = '../../../out/FM_Burner/'

fuel_name = ['C2H4', 'C3H6', 'C3H8', 'CH4']
res_name = ['2cm', '1cm', '5mm']

# --- First part: XO2, eta, Chi_R ---
for i_fuel in range(4):
    for fds_resolution in range(3):
        dev_file = os.path.join(outdir, f"FM_15cm_Burner_{fuel_name[i_fuel]}_{res_name[fds_resolution]}_devc.csv")
        hrr_file = os.path.join(outdir, f"FM_15cm_Burner_{fuel_name[i_fuel]}_{res_name[fds_resolution]}_hrr.csv")

        DEV = pd.read_csv(dev_file, header=1)
        HRR = pd.read_csv(hrr_file, header=1)

        Time_FDS = DEV['Time'].to_numpy()
        XO2_FDS = DEV['XO2'].to_numpy()
        Qdot_FDS = HRR['HRR'].to_numpy()
        Qrad_FDS = HRR['Q_RADI'].to_numpy()
        ntp = len(Time_FDS)

        out_file = os.path.join(outdir, f"FM_15cm_Burner_{fuel_name[i_fuel]}_{res_name[fds_resolution]}.csv")
        with open(out_file, 'w') as fid:
            fid.write('XO2,eta,Chi_R\n')
            for ii in range(ntp):
                eta = Qdot_FDS[ii] / np.max(Qdot_FDS)
                chi_r = max(0, -Qrad_FDS[ii] / max(0.001, Qdot_FDS[ii]))
                fid.write(f"{XO2_FDS[ii]:5.3f},{eta:6.2f},{chi_r:6.2f}\n")

# --- Second part: Chi_r vs Time for different O2 levels ---
O2_name = ['20p9', '19p0', '16p8', '15p2']

for i_O2 in range(4):
    for fds_resolution in range(3):
        hrr_file = os.path.join(outdir, f"FM_15cm_Burner_C2H4_{O2_name[i_O2]}_{res_name[fds_resolution]}_hrr.csv")
        HRR = pd.read_csv(hrr_file, header=1)

        Time_FDS = HRR['Time'].to_numpy()
        Qdot_FDS = HRR['HRR'].to_numpy()
        Qrad_FDS = HRR['Q_RADI'].to_numpy()
        ntp = len(Time_FDS)

        chi_r = np.maximum(0, -Qrad_FDS / np.maximum(0.001, Qdot_FDS))

        out_file = os.path.join(outdir, f"FM_15cm_Burner_C2H4_{O2_name[i_O2]}_{res_name[fds_resolution]}_chir.csv")
        with open(out_file, 'w') as fid:
            fid.write('Time,Chi_r\n')
            for ii in range(1, ntp):  # start from 2nd index like MATLAB
                fid.write(f"{Time_FDS[ii]:5.3f},{chi_r[ii]:6.3f}\n")
