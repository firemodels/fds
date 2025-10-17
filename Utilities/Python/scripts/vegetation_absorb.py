
# Convert last row of devc file for WUI/vegetation_absorb into a column to be read by dataplot

import pandas as pd

outdir = '../../Verification/WUI/'

M = pd.read_csv(outdir + 'vegetation_absorb_devc.csv', skiprows=2)

with open(outdir + 'vegetation_absorb_FDS.csv', 'w', newline='') as fid:
    fid.write('mpuv,rad\n')
    mpuv = [0.0, 0.1, 0.2, 0.4, 0.8, 1.6]
    for i in range(6):
        fid.write(f'{mpuv[i]:5.1f},{M.iloc[-1, i+1]:5.2f}\n')

