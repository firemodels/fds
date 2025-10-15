# Roger Wang
# 7-20-11
# ashrae_7.m
#
# Converted by Floyd

import os
import numpy as np
import pandas as pd

datadir =  '../../Verification/HVAC/'
filename = [
    'ashrae_7_exp.csv',
    'ashrae7_fixed_flow_devc.csv',
    'ashrae7_quadratic_devc.csv',
    'ashrae7_table_devc.csv'
]

label = [
    'Experiment &',
    'Fixed Flow & ',
    'Quadratic & ',
    'Table & '
]

# duct = ['1', '2', '3', '4', '5', '56', '6', '7']

pressure = np.zeros((4, 8))

skip_case = False

for i in range(len(filename)):
   name = datadir+filename[i]
   if not os.path.exists(name):
      skip_case = True
      print('Error: File ', filename[i], ' does not exist. Skipping case.')

if skip_case: quit()

M = pd.read_csv(datadir+filename[0],skiprows=6,header=None)

for n in range(1, len(filename)):
   df_sim = pd.read_csv(datadir+filename[n],skiprows=2,header=None)
   last_row_data = df_sim.iloc[-1].values
   pressure[n, :] = last_row_data[1:9]

# --- Write to LaTeX File ---

texname = '../../Manuals/FDS_Verification_Guide/SCRIPT_FIGURES/ashrae_7.tex'

with open(texname, 'w') as fid:

   # Boilerplate LaTeX environment opening
   fid.write('\\begin{center}\n')
   fid.write('\\begin{tabular}{|c|c|c|c|c|c|c|c|c|} \\hline\n')

   # Header Row
   fid.write('Duct Number & 1 & 2 & 3 & 4 & 5 & 56 & 6 & 7 \\\\ \\hline\n')

   # Experiment Row (M)
   format_str_exp =label[0] + ' %6.3f &' * 7 + ' %6.3f \\\\ \n'

   # The '*' operator unpacks the NumPy array 'M' into arguments for the format string.
   fid.write(format_str_exp % tuple(M.iloc[0,1:9]))

   # Simulation Rows (pressure)
   for i in range(1,  len(filename)):

      # Build the format string using the label prefix
      format_str_sim = label[i] + ('%6.3f & ' * 7) + '%6.3f \\\\ \n'

      # Unpack the pressure data and write the line
      fid.write(format_str_sim % tuple(pressure[i,:]))

   # Boilerplate LaTeX environment closing
   fid.write('\\hline\n')
   fid.write('\\end{tabular}\n')
   fid.write('\\end{center}\n')

print(f"\n? Successfully generated LaTeX table in {texname}")