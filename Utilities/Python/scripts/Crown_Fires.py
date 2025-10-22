"""
McGrattan
10-29-2019
Crown_Fires.py

Read the Crown_Fires *_cat_devc.csv files and determine the rate of spread based on the time history of front position.
Write the results to a file that will be plotted via dataplot.m.
"""

import numpy as np
import pandas as pd
import os

# Parameters
outdir = '../../../out/Crown_Fires/'

# Get only *_cat_devc.csv files (like MATLAB dir)
file_list = [f for f in os.listdir(outdir) if f.endswith('_cat_devc.csv')]
file_list.sort()  # optional, to keep consistent ordering with MATLAB

wind_speed = []
slope = []

for fname in file_list:
    full_path = os.path.join(outdir, fname)

    # Read CSV (skip 2 header rows like MATLAB importdata)
    M = pd.read_csv(full_path, skiprows=2)
    M_data = M.to_numpy()

    # Extract rows satisfying conditions (700<=col2<=900, 30<col1<300)
    indices = np.where(
        (M_data[:, 1] >= 700) &
        (M_data[:, 1] <= 900) &
        (M_data[:, 0] > 30) &
        (M_data[:, 0] < 300)
    )[0]

    # Mean wind speed (col3)
    wind_speed.append(np.mean(M_data[indices, 2]))

    # Polyfit slope for col2 vs col1
    p = np.polyfit(M_data[indices, 0], M_data[indices, 1], 1)
    slope.append(p[0])

# Write output file
with open(os.path.join(outdir, 'ROS.csv'), 'w') as fid:
    fid.write('km/h,m/min\n')
    fid.write('U,ROS\n')
    for u, s in zip(wind_speed, slope):
        fid.write(f'{3.6*u:4.1f},{60*s:6.2f}\n')
