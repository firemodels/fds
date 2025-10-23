
# FDS Verification Guide, write convergence info to Verification/Heat_Transfer/convective_cooling_error.csv

import pandas as pd
import os

outdir = '../../Verification/Heat_Transfer/'

infile = [outdir + 'convective_cooling_p1_devc.csv',
          outdir + 'convective_cooling_p05_devc.csv',
          outdir + 'convective_cooling_p025_devc.csv',
          outdir + 'convective_cooling_p01_devc.csv',
          outdir + 'convective_cooling_p005_devc.csv',
          outdir + 'convective_cooling_p0025_devc.csv',
          outdir + 'convective_cooling_p00125_devc.csv']

dx = [0.1, 0.05, 0.025, 0.01, 0.005, 0.0025, 0.00125]
T = [0.0] * 7

for k in range(7):
    if not os.path.exists(infile[k]):
        print(f'Error: File {infile[k]} does not exist. Skipping case.')
        exit()

    M_10 = pd.read_csv(infile[k], skiprows=2, header=None)
    T[k] = M_10.iloc[-1, 1]

with open(outdir + 'convective_cooling_error.csv', 'w') as fid11:
    fid11.write('dx,dx^2,relative error\n')
    T_exact = 295.3011488157
    for j in range(7):
        relative_error = abs(T[j] - T_exact) / T_exact
        fid11.write(f'{0.01*dx[j]:8.3e}, {0.1*dx[j]**2:8.3e}, {relative_error:9.5e}\n')

