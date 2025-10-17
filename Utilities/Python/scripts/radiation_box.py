# radiation_box.py     Collect heat fluxes from radiation box cases to a single csv file in Verification/Radiation folder

import os
import numpy as np
import pandas as pd

# Define positions
x = np.array([
    2.50E-02, 7.50E-02, 1.25E-01, 1.75E-01, 2.25E-01,
    2.75E-01, 3.25E-01, 3.75E-01, 4.25E-01, 4.75E-01,
    5.25E-01, 5.75E-01, 6.25E-01, 6.75E-01, 7.25E-01,
    7.75E-01, 8.25E-01, 8.75E-01, 9.25E-01, 9.75E-01
])

# Directory
dir = '../../Verification/Radiation/'

# Input files
infile = [
    f'{dir}radiation_box__20___50_devc.csv',
    f'{dir}radiation_box__20__100_devc.csv',
    f'{dir}radiation_box__20__300_devc.csv',
    f'{dir}radiation_box__20_1000_devc.csv',
    f'{dir}radiation_box__20_2000_devc.csv',
    f'{dir}radiation_box_100___50_devc.csv',
    f'{dir}radiation_box_100__100_devc.csv',
    f'{dir}radiation_box_100__300_devc.csv',
    f'{dir}radiation_box_100_1000_devc.csv',
    f'{dir}radiation_box_100_2000_devc.csv'
]

# Corresponding labels
label = [
    'Flux_20_50', 'Flux_20_100', 'Flux_20_300', 'Flux_20_1000', 'Flux_20_2000',
    'Flux_100_50', 'Flux_100_100', 'Flux_100_300', 'Flux_100_1000', 'Flux_100_2000'
]

# Initialize matrix for flux data
flux = np.zeros((20, 10))

# Read data and populate flux
for n in range(10):
    if not os.path.exists(infile[n]):
        print(f'Error: File {infile[n]} does not exist. Skipping case.')
    else:
        # Read CSV file starting from the 4th row
        data = pd.read_csv(infile[n], skiprows=3, header=None)
        t = data.iloc[0, 0]
        flux[:, n] = data.iloc[0, 1:21].values

# Output file
filename = f'{dir}radiation_box_devc.csv'

# Write to CSV
with open(filename, 'w') as fid:
    fid.write('Position, ' + ', '.join(label) + '\n')
    for i in range(20):
        fid.write(f'{x[i]:.6f}, ' + ', '.join(f'{flux[i, j]:.6f}' for j in range(10)) + '\n')
