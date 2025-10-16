# C Weinschenk
# flame_species.m
# 10-2011
# Combine outputs from Methane_flame_lumped and
# Methane_flame_primitive into 1 file
# ------------------------------------------------
#
# Converted by Floyd
# 10/16/2025

import pandas as pd
import numpy as np
import os
import csv
from typing import List, Tuple, Optional


# Write directory
outdir = '../../Verification/Species/'
filename = [
    'methane_flame_primitive_devc.csv',
    'methane_flame_primitive_2_devc.csv',
    'methane_flame_lumped_devc.csv',
    'methane_flame_lumped_fuel_devc.csv',
    'methane_flame_lumped_ox_devc.csv',
]

skip_case = False

for i in range(len(filename)):
   name = outdir+filename[i]
   if not os.path.exists(name):
      skip_case = True
      print('Error: File ', filename[i], ' does not exist. Skipping case.')

if skip_case: quit()

h_mf_p = pd.read_csv(outdir+filename[0],skiprows=1,nrows=1,header=None)
h_mf_p_2 = pd.read_csv(outdir+filename[1],skiprows=1,nrows=1,header=None)
h_mf_l = pd.read_csv(outdir+filename[2],skiprows=1,nrows=1,header=None)
h_mf_l_f = pd.read_csv(outdir+filename[3],skiprows=1,nrows=1,header=None)
h_mf_l_o = pd.read_csv(outdir+filename[4],skiprows=1,nrows=1,header=None)

mf_p = pd.read_csv(outdir+filename[0],skiprows=2,header=None)
mf_p_2 = pd.read_csv(outdir+filename[1],skiprows=2,header=None)
mf_l = pd.read_csv(outdir+filename[2],skiprows=2,header=None)
mf_l_f = pd.read_csv(outdir+filename[3],skiprows=2,header=None)
mf_l_o = pd.read_csv(outdir+filename[4],skiprows=2,header=None)

header_lp = tuple(h_mf_p.iloc[0,:]) + tuple(h_mf_l.iloc[0,1:])

print(header_lp)

data_combined_lp = pd.concat([mf_p, mf_l.iloc[:, 1:4]], axis=1)
data_combined_lp.columns = header_lp

print(data_combined_lp)

header_ml =  tuple(h_mf_p_2.iloc[0,:]) +  tuple(h_mf_l_f.iloc[0,1:3]) +  tuple(h_mf_l_o.iloc[0,1:3])
print(header_ml)
data_combined_ml = pd.concat([mf_p_2, mf_l_f.iloc[:, 1:3], mf_l_o.iloc[:, 1:3]], axis=1)
data_combined_ml.columns = header_ml

print(data_combined_ml)

output_file_1 = outdir+'methane_flame_lumpedprimitive2.csv'

# The MATLAB script uses a hardcoded units line for the final output.
hardcoded_units = ['s', 'kg', 'kg', 'kg', 'kg', 'kg', 'kg']

# Create a temporary DataFrame for the two header rows
header_df_1 = pd.DataFrame([hardcoded_units, header_lp], columns=data_combined_lp.columns)

# Concatenate headers and data
df_out_1 = pd.concat([header_df_1, data_combined_lp]).reset_index(drop=True)

# Write to CSV without the pandas index or automatically generated header
df_out_1.to_csv(output_file_1, header=False, index=False)
print(f"Successfully wrote combined data to {output_file_1}")

output_file_2 = outdir+'methane_flame_multilumped2.csv'

# Create a temporary DataFrame for the two header rows
header_df_2 = pd.DataFrame([hardcoded_units, header_ml], columns=data_combined_ml.columns)

# Concatenate headers and data
df_out_2 = pd.concat([header_df_2, data_combined_ml]).reset_index(drop=True)

# Write to CSV without the pandas index or automatically generated header
df_out_2.to_csv(output_file_2, header=False, index=False)
print(f"Successfully wrote multi-lumped data to {output_file_2}")
