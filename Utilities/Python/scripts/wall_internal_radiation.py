# wall_internal_radiation.py

import os
import pandas as pd

# File paths
infile = '../../Verification/Radiation/wall_internal_radiation_devc.csv'
outfile = '../../Manuals/FDS_Verification_Guide/SCRIPT_FIGURES/wall_internal_radiation.tex'

# Check if the input file exists
if not os.path.exists(infile):
    print(f'Error: File {infile} does not exist. Skipping case.')
else:
    # Read CSV file starting from the 4th row
    data = pd.read_csv(infile, skiprows=3, header=None)
    t = data.iloc[0, 0]
    flux = data.iloc[0, 1:6].values
    
    # Open the output file and write LaTeX content
    with open(outfile, 'w') as fid:
        fid.write('\\begin{center}\n')
        fid.write('\\begin{tabular}{|c|c|c|} \\hline\n')
        fid.write('$\\tau$      & $S(\\tau)$   & FDS \\\\\n')
        fid.write('            & (kW/m$^2$)  & (kW/m$^2$) \\\\ \\hline\\hline\n')
        fid.write(f'0.01        & 2.897       & {-flux[0]:6.3f} \\\\\n')
        fid.write(f'0.1         & 24.94       & {-flux[1]:5.2f} \\\\\n')
        fid.write(f'0.5         & 82.95       & {-flux[2]:5.2f} \\\\\n')
        fid.write(f'1.0         & 116.3       & {-flux[3]:5.1f} \\\\\n')
        fid.write(f'10.         & 149.0       & {-flux[4]:5.1f} \\\\ \\hline\n')
        fid.write('\\end{tabular}\n')
        fid.write('\\end{center}\n')

