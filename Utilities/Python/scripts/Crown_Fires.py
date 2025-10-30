
# Read the Crown_Fires *_cat_devc.csv files and determine the rate of spread based on the time history of front position.
# Write the results to a file that will be plotted via dataplot.py.

import numpy as np
import pandas as pd
import os

outdir = '../../../out/Crown_Fires/'

file_list = [f for f in os.listdir(outdir) if f.endswith('_cat_devc.csv')]
file_list.sort() 

wind_speed = []
slope = []

for fname in file_list:
    full_path = os.path.join(outdir, fname)

    M = pd.read_csv(full_path, skiprows=2, header=None)
    M_data = M.to_numpy()

    # Extract rows satisfying conditions (700<=x<=900, 30<Time<300)
    indices = np.where(
        (M_data[:, 1] >= 700) &
        (M_data[:, 1] <= 900) &
        (M_data[:, 0] > 30) &
        (M_data[:, 0] < 300)
    )[0]

    # Mean wind speed (U10)
    wind_speed.append(np.mean(M_data[indices, 2]))

    # Polyfit slope for x vs Time
    p = np.polyfit(M_data[indices, 0], M_data[indices, 1], 1)
    slope.append(p[0])

with open(os.path.join(outdir, 'ROS.csv'), 'w') as fid:
    fid.write('km/h,m/min\n')
    fid.write('U,ROS\n')
    for u, s in zip(wind_speed, slope):
        fid.write(f'{3.6*u:4.1f},{60*s:6.2f}\n')

