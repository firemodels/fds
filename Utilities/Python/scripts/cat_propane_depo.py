# Overholt
# 6-13-2012
# cat_propane_depo.m
#
# Concatenates columns from Propane flame deposition FDS cases (/Verification/Aerosols)
#
# Converted by Floyd
# 10-16-2025


import pandas as pd
import numpy as np
import os


outdir = '../../Verification/Aerosols/'

# --- Condensed phase aerosol (Wall Deposition) ---

# List of files for condensed phase aerosol (wall deposition)
wall_files = [
    'propane_flame_deposition_devc.csv',
    'propane_flame_deposition_none_devc.csv',
    'propane_flame_deposition_gravitational_devc.csv',
    'propane_flame_deposition_thermophoretic_devc.csv',
    'propane_flame_deposition_turbulent_devc.csv'
]

wall_columns = [
    'depo_all',
    'depo_none',
    'depo_gravitational',
    'depo_thermophoretic',
    'depo_turbulent'
]

skip_case = False

for i in range(len(wall_files)):
   name = outdir+wall_files[i]
   if not os.path.exists(name):
      skip_case = True
      print('Error: File ', wall_files[i], ' does not exist. Skipping case.')

if skip_case: quit()

# Read data from the wall deposition files
wall_data = {}
for i, file in enumerate(wall_files):
   wall_data[wall_columns[i]] = pd.read_csv(outdir+wall_files[i], skiprows=2, header=None)

# Extract the relevant columns:
time_col = wall_data[wall_columns[0]].iloc[:, 0]
depo_all_col = wall_data[wall_columns[0]].iloc[:, 1]
depo_none_col = wall_data[wall_columns[1]].iloc[:, 1]
depo_grav_col = wall_data[wall_columns[2]].iloc[:, 1]
depo_therm_col = wall_data[wall_columns[3]].iloc[:, 1]
depo_turb_col = wall_data[wall_columns[4]].iloc[:, 1]

# Create a DataFrame D1
D1 = pd.DataFrame({
    'Time': time_col,
    'depo_all': depo_all_col,
    'depo_none': depo_none_col,
    'depo_gravitational': depo_grav_col,
    'depo_thermophoretic': depo_therm_col,
    'depo_turbulent': depo_turb_col
})

# Add the two header rows defined in MATLAB
H1_units = ['s', 'kg', 'kg', 'kg', 'kg', 'kg']
H1_names = ['Time', 'depo_all', 'depo_none', 'depo_gravitational', 'depo_thermophoretic', 'depo_turbulent']

# Create a dummy DataFrame for the headers
H1_df = pd.DataFrame([H1_units, H1_names], columns=D1.columns)

# Concatenate headers and data
D1_final = pd.concat([H1_df, D1]).reset_index(drop=True)

D1_final.to_csv(outdir+'propane_flame_deposition_cat_wall.csv', header=False, index=False)
print(f"Written condensed phase aerosol data to propane_flame_deposition_cat_wall.csv")

# List of files for gas phase aerosol (mass loss)
gas_files = [
    'propane_flame_deposition_mass.csv',
    'propane_flame_deposition_none_mass.csv',
    'propane_flame_deposition_gravitational_mass.csv',
    'propane_flame_deposition_thermophoretic_mass.csv',
    'propane_flame_deposition_turbulent_mass.csv'
]

gas_columns = [
    'depo_all', # Reused names for D2 columns, N1, N2, etc. are the source files
    'depo_none',
    'depo_gravitational',
    'depo_thermophoretic',
    'depo_turbulent'
]

for i in range(len(gas_files)):
   name = outdir+gas_files[i]
   if not os.path.exists(name):
      skip_case = True
      print('Error: File ', gas_files[i], ' does not exist. Skipping case.')

if skip_case: quit()

# Read data from the gas mass files
gas_data = {}
for i, file in enumerate(gas_files):
   gas_data[gas_columns[i]] = pd.read_csv(outdir+gas_files[i], skiprows=2, header=None)

# Extract the relevant columns:
N1_time_col = gas_data[gas_columns[0]].iloc[:, 0]
N1_mass_col = gas_data[gas_columns[0]].iloc[:, 7]
N2_mass_col = gas_data[gas_columns[1]].iloc[:, 7]
N3_mass_col = gas_data[gas_columns[2]].iloc[:, 7]
N4_mass_col = gas_data[gas_columns[3]].iloc[:, 7]
N5_mass_col = gas_data[gas_columns[4]].iloc[:, 7]

# Create a DataFrame D2
D2 = pd.DataFrame({
    'Time': N1_time_col,
    'depo_all': N1_mass_col,
    'depo_none': N2_mass_col,
    'depo_gravitational': N3_mass_col,
    'depo_thermophoretic': N4_mass_col,
    'depo_turbulent': N5_mass_col
})

# Add the two header rows defined in MATLAB (H2 is identical to H1)
H2_units = ['s', 'kg', 'kg', 'kg', 'kg', 'kg']
H2_names = ['Time', 'depo_all', 'depo_none', 'depo_gravitational', 'depo_thermophoretic', 'depo_turbulent']

# Create a dummy DataFrame for the headers
H2_df = pd.DataFrame([H2_units, H2_names], columns=D2.columns)

# Concatenate headers and data
D2_final = pd.concat([H2_df, D2]).reset_index(drop=True)

# Write to CSV
D2_final.to_csv(outdir+'propane_flame_deposition_cat_gas.csv', header=False, index=False)
print(f"Written gas phase aerosol data to propane_flame_deposition_cat_gas.csv")

# Create DataFrame D3
D3 = pd.DataFrame({
    'Time': N1_time_col,
    'depo_all': depo_all_col + N1_mass_col,
    'depo_none': depo_none_col + N2_mass_col,
    'depo_gravitational':depo_grav_col + N3_mass_col,
    'depo_thermophoretic': depo_therm_col + N4_mass_col,
    'depo_turbulent': depo_turb_col + N5_mass_col
})

# Add the two header rows defined in MATLAB (H3 is identical to H1 and H2)
H3_units = ['s', 'kg', 'kg', 'kg', 'kg', 'kg']
H3_names = ['Time', 'depo_all', 'depo_none', 'depo_gravitational', 'depo_thermophoretic', 'depo_turbulent']

# Create a dummy DataFrame for the headers
H3_df = pd.DataFrame([H3_units, H3_names], columns=D3.columns)

# Concatenate headers and data
D3_final = pd.concat([H3_df, D3]).reset_index(drop=True)

# Write to CSV
D3_final.to_csv(outdir+ 'propane_flame_deposition_cat_total.csv', header=False, index=False)
print(f"Written total aerosol data to propane_flame_deposition_cat_total.csv")