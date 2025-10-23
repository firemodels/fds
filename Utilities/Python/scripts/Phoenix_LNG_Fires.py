
# Reformat heat flux output for Phoenix_LNG_Fires

import os
import numpy as np
import pandas as pd

outdir = '../../../out/Phoenix_LNG_Fires/'

def write_csv(filename, header, rows):
    """Write a CSV with header and rows, using %5.1f formatting."""
    with open(filename, 'w') as f:
        f.write(header + '\n')
        for row in rows:
            f.write(','.join(f"{val:5.1f}" for val in row) + '\n')

# --- Phoenix01 ---
DEV = pd.read_csv(os.path.join(outdir, 'Phoenix01_devc.csv'), skiprows=2, header=None)

rows_hf1 = [
    [128.6, DEV.iloc[12,1], 108.2, DEV.iloc[12,4],  92.0, DEV.iloc[12,7], 112.4, DEV.iloc[12,10]],
    [178.5, DEV.iloc[12,2], 157.4, DEV.iloc[12,5], 141.2, DEV.iloc[12,8], 162.3, DEV.iloc[12,11]],
    [228.3, DEV.iloc[12,3], 207.3, DEV.iloc[12,6], 191.1, DEV.iloc[12,9], 212.1, DEV.iloc[12,12]],
]
write_csv(os.path.join(outdir, 'Phoenix01_HF.csv'),
          'Pos-North-Wide,North-Wide,Pos-East-Wide,East-Wide,Pos-South-Wide,South-Wide,Pos-West-Wide,West-Wide',
          rows_hf1)

rows_nahf1 = [
    [  7.6, DEV.iloc[12,13],  7.4, DEV.iloc[12,20]],
    [ 21.2, DEV.iloc[12,14], 23.1, DEV.iloc[12,21]],
    [ 35.4, DEV.iloc[12,15], 35.3, DEV.iloc[12,22]],
    [ 49.2, DEV.iloc[12,16], 49.2, DEV.iloc[12,23]],
    [ 62.9, DEV.iloc[12,17], 65.0, DEV.iloc[12,24]],
    [174.4, DEV.iloc[12,18], 39.4, DEV.iloc[12,25]],
    [134.0, DEV.iloc[12,19], 36.3, DEV.iloc[12,26]],
]
write_csv(os.path.join(outdir, 'Phoenix01_NAHF.csv'),
          'North_Height,North_Flux,South_Height,South_Flux',
          rows_nahf1)

# --- Phoenix02 ---
DEV = pd.read_csv(os.path.join(outdir, 'Phoenix02_devc.csv'), skiprows=2, header=None)

rows_hf2 = [
    [133.0, DEV.iloc[12,1], 103.4, DEV.iloc[12,4],  87.6, DEV.iloc[12,7], 117.7, DEV.iloc[12,10]],
    [182.9, DEV.iloc[12,2], 152.0, DEV.iloc[12,5], 136.8, DEV.iloc[12,8], 167.5, DEV.iloc[12,11]],
    [232.7, DEV.iloc[12,3], 202.1, DEV.iloc[12,6], 186.7, DEV.iloc[12,9], 217.5, DEV.iloc[12,12]],
]
write_csv(os.path.join(outdir, 'Phoenix02_HF.csv'),
          'Pos-North-Wide,North-Wide,Pos-East-Wide,East-Wide,Pos-South-Wide,South-Wide,Pos-West-Wide,West-Wide',
          rows_hf2)

rows_nahf2 = [
    [ 15.0, DEV.iloc[12,13],  15.0, DEV.iloc[12,20]],
    [ 30.0, DEV.iloc[12,14],  30.0, DEV.iloc[12,21]],
    [ 55.0, DEV.iloc[12,15],  55.0, DEV.iloc[12,22]],
    [ 85.0, DEV.iloc[12,16],  85.0, DEV.iloc[12,23]],
    [120.0, DEV.iloc[12,17], 120.0, DEV.iloc[12,24]],
    [120.0, DEV.iloc[12,18], 120.0, DEV.iloc[12,25]],
    [120.0, DEV.iloc[12,19], 120.0, DEV.iloc[12,26]],
]
write_csv(os.path.join(outdir, 'Phoenix02_NAHF.csv'),
          'North_Height,North_Flux,South_Height,South_Flux',
          rows_nahf2)

