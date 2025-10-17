
# This script must be run prior to dataplot both to generate the expected results
# and to normalize the PDPA results from the FDS line file.

import numpy as np
import pandas as pd
import os

outdir = '../../Verification/WUI/'

z0 = 5
nz = 80
L = 10
dz = L/nz
z = np.linspace(dz/2, L-dz/2, nz)
zp = z - z0
RHO = 1.199  # kg/m3
D_1 = 0.1/RHO  # 0.0834 m2/s
D_2 = 0.01/RHO  # 0.00834 m2/s
t = 15  # s

f_1 = 1/(4*np.pi*D_1*np.sqrt(t))*np.exp(-zp**2/(4*D_1*t))
f_1 = f_1/np.sum(f_1*dz)  # normalize

f_2 = 1/(4*np.pi*D_2*np.sqrt(t))*np.exp(-zp**2/(4*D_2*t))
f_2 = f_2/np.sum(f_2*dz)  # normalize

# print expected results to data directory

fid = open(outdir + 'random_walk.csv', 'wt')
fid.write('{}, {}, {}\n'.format('z', 'f_1', 'f_2'))
for i in range(len(z)):
    fid.write('{:f}, {:f}, {:f}\n'.format(z[i], f_1[i], f_2[i]))
fid.close()

# pdpa_radius = 0.5; % m
# pdpa_volume = 4/3*pi*pdpa_radius^3; % m^3

# normalize and overwrite the FDS results for random_walk_1

filename = outdir + 'random_walk_1_line.csv'

if not os.path.exists(filename):  # skip_case_if

    print('Error: File ' + filename + ' does not exist. Skipping case.')

else:

    M = pd.read_csv(outdir + 'random_walk_1_line.csv', skiprows=1)
    npv_z = M['npv-z'].values
    dz = (np.max(npv_z) - np.min(npv_z))/(len(npv_z) - 1)
    npv = M['npv'].values
    nt = np.sum(npv)

    fid = open(outdir + 'random_walk_1.csv', 'wt')
    fid.write('{}, {}\n'.format('m', ' '))
    fid.write('{}, {}\n'.format('npv-z', 'npv'))
    for i in range(len(npv_z)):
        fid.write('{:f}, {:f}\n'.format(npv_z[i], npv[i]/nt/dz))
    fid.close()

# normalize and overwrite the FDS results for random_walk_2

filename = outdir + 'random_walk_2_line.csv'

if not os.path.exists(filename):  # skip_case_if

    print('Error: File ' + filename + ' does not exist. Skipping case.')

else:

    M = pd.read_csv(outdir + 'random_walk_2_line.csv', skiprows=1)
    npv_z = M['npv-z'].values
    dz = (np.max(npv_z) - np.min(npv_z))/(len(npv_z) - 1)
    npv = M['npv'].values
    nt = np.sum(npv)

    fid = open(outdir + 'random_walk_2.csv', 'wt')
    fid.write('{}, {}\n'.format('m', ' '))
    fid.write('{}, {}\n'.format('npv-z', 'npv'))
    for i in range(len(npv_z)):
        fid.write('{:f}, {:f}\n'.format(npv_z[i], npv[i]/nt/dz))
    fid.close()

