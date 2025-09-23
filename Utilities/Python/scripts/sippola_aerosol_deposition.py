#!/usr/bin/env python
"""
Converted from

Overholt
10-30-2012
sippola_aerosol_deposition.m

Calculations for Aerosol Deposition FDS Validation cases
located in (/Validation/Sippola_Aerosol_Deposition)
"""

import os
import numpy as np

# Output directory
outdir = '../../../out/Sippola_Aerosol_Deposition/'

# Ensure output directory exists
if not os.path.exists(outdir):
    os.makedirs(outdir)

# Filenames
filenames = [f'Sippola_Test_{i:02d}_devc.csv' for i in range(1, 32)]

# Friction velocity (m/s)
friction_velocity = [
    0.12, 0.12, 0.12, 0.13, 0.12, 0.28, 0.26, 0.26,
    0.27, 0.28, 0.28, 0.45, 0.42, 0.44, 0.46, 0.45,
    0.16, 0.16, 0.16, 0.16, 0.16, 0.37, 0.37, 0.37,
    0.38, 0.38, 0.62, 0.62, 0.62, 0.64, 0.64
]

# Air velocity (m/s)
air_velocity = [
    2.2, 2.2, 2.1, 2.2, 2.2, 5.3, 5.2, 5.2,
    5.4, 5.3, 5.3, 9.0, 9.0, 8.8, 9.2, 9.1,
    2.2, 2.2, 2.2, 2.2, 2.2, 5.3, 5.2, 5.2,
    5.3, 5.3, 8.9, 8.7, 8.9, 8.9, 8.9
]

# Particle diameters (µm)
particle_sizes = [1, 3, 5, 9, 16]

# Arrays for deposition
Vd_ceiling = np.zeros((len(filenames), len(particle_sizes)))
Vd_wall    = np.zeros((len(filenames), len(particle_sizes)))
Vd_floor   = np.zeros((len(filenames), len(particle_sizes)))

# Read CSVs and extract deposition
for i, fname in enumerate(filenames):
    file_path = os.path.join(outdir, fname)
    if not os.path.exists(file_path):
        continue  # skip missing files silently
    
    data = np.genfromtxt(file_path, delimiter=',', skip_header=2)
    
    for j, size in enumerate(particle_sizes):
        Vd_ceiling[i, j] = data[j, 6]
        Vd_wall[i, j]    = data[j, 7]
        Vd_floor[i, j]   = data[j, 5]

# Non-dimensional deposition
Vd_ceiling_nd = Vd_ceiling / np.array(friction_velocity)[:, None]
Vd_wall_nd    = Vd_wall / np.array(friction_velocity)[:, None]
Vd_floor_nd   = Vd_floor / np.array(friction_velocity)[:, None]

# Write output table
output_file = os.path.join(outdir, 'Sippola_All_Tests.csv')
with open(output_file, 'w') as fid:
    header = ['Test', 'AirVelocity', 'FrictionVelocity']
    for size in particle_sizes:
        header += [f'Ceiling_{size}um', f'Wall_{size}um', f'Floor_{size}um']
    for size in particle_sizes:
        header += [f'Ceiling_{size}um_ND', f'Wall_{size}um_ND', f'Floor_{size}um_ND']
    
    fid.write(','.join(header) + '\n')
    
    for i in range(len(filenames)):
        row = [str(i+1), f"{air_velocity[i]:.2f}", f"{friction_velocity[i]:.2f}"]
        for j in range(len(particle_sizes)):
            row += [f"{Vd_ceiling[i,j]:.3e}", f"{Vd_wall[i,j]:.3e}", f"{Vd_floor[i,j]:.3e}"]
        for j in range(len(particle_sizes)):
            row += [f"{Vd_ceiling_nd[i,j]:.3e}", f"{Vd_wall_nd[i,j]:.3e}", f"{Vd_floor_nd[i,j]:.3e}"]
        
        fid.write(','.join(row) + '\n')
