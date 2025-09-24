#!/usr/bin/env python3
"""
Jason Floyd
Ranz_Marshall.py
"""

import os
import pandas as pd

outdir = '../../../out/Ranz_Marshall/'

def process_table(filenames, output_file):
    """Check files, extract last row col6, write results."""
    # Check for missing files
    missing = [f for f in filenames if not os.path.exists(os.path.join(outdir, f))]
    if missing:
        for f in missing:
            print(f"Error: File {os.path.join(outdir, f)} does not exist. Skipping case.")
        return  # same as MATLAB return

    d2dt = []
    for f in filenames:
        df = pd.read_csv(os.path.join(outdir, f), skiprows=2)
        d2dt.append(df.iloc[-1, 5])  # last row, column index 5 (6th column)

    # Write out results file
    with open(os.path.join(outdir, output_file), 'w') as fid:
        fid.write('Test,dD2dt\n')
        for i, val in enumerate(d2dt, start=1):
            fid.write(f'{i},  {val:3.1e}\n')

# Table 1
filenames1 = [
    'Ranz_Marshall_Table_1_1_devc.csv', 'Ranz_Marshall_Table_1_2_devc.csv',
    'Ranz_Marshall_Table_1_3_devc.csv', 'Ranz_Marshall_Table_1_4_devc.csv',
    'Ranz_Marshall_Table_1_5_devc.csv', 'Ranz_Marshall_Table_1_6_devc.csv',
    'Ranz_Marshall_Table_1_7_devc.csv', 'Ranz_Marshall_Table_1_8_devc.csv',
    'Ranz_Marshall_Table_1_9_devc.csv', 'Ranz_Marshall_Table_1_10_devc.csv',
    'Ranz_Marshall_Table_1_11_devc.csv', 'Ranz_Marshall_Table_1_12_devc.csv',
    'Ranz_Marshall_Table_1_13_devc.csv', 'Ranz_Marshall_Table_1_14_devc.csv',
    'Ranz_Marshall_Table_1_15_devc.csv', 'Ranz_Marshall_Table_1_16_devc.csv',
    'Ranz_Marshall_Table_1_17_devc.csv', 'Ranz_Marshall_Table_1_18_devc.csv',
    'Ranz_Marshall_Table_1_19_devc.csv'
]
process_table(filenames1, 'Ranz_Marshall_Table_1.csv')

# Table 2
filenames2 = [
    'Ranz_Marshall_Table_2_1_devc.csv', 'Ranz_Marshall_Table_2_2_devc.csv',
    'Ranz_Marshall_Table_2_3_devc.csv', 'Ranz_Marshall_Table_2_4_devc.csv',
    'Ranz_Marshall_Table_2_5_devc.csv', 'Ranz_Marshall_Table_2_6_devc.csv',
    'Ranz_Marshall_Table_2_7_devc.csv', 'Ranz_Marshall_Table_2_8_devc.csv',
    'Ranz_Marshall_Table_2_9_devc.csv'
]
process_table(filenames2, 'Ranz_Marshall_Table_2.csv')

# Table 3
filenames3 = [
    'Ranz_Marshall_Table_3_1_devc.csv', 'Ranz_Marshall_Table_3_2_devc.csv',
    'Ranz_Marshall_Table_3_3_devc.csv', 'Ranz_Marshall_Table_3_4_devc.csv',
    'Ranz_Marshall_Table_3_5_devc.csv', 'Ranz_Marshall_Table_3_6_devc.csv',
    'Ranz_Marshall_Table_3_7_devc.csv', 'Ranz_Marshall_Table_3_8_devc.csv',
    'Ranz_Marshall_Table_3_9_devc.csv'
]
process_table(filenames3, 'Ranz_Marshall_Table_3.csv')

# Table 4
filenames4 = [
    'Ranz_Marshall_Table_4_1_devc.csv', 'Ranz_Marshall_Table_4_2_devc.csv',
    'Ranz_Marshall_Table_4_3_devc.csv', 'Ranz_Marshall_Table_4_4_devc.csv',
    'Ranz_Marshall_Table_4_5_devc.csv', 'Ranz_Marshall_Table_4_6_devc.csv',
    'Ranz_Marshall_Table_4_7_devc.csv', 'Ranz_Marshall_Table_4_8_devc.csv',
    'Ranz_Marshall_Table_4_9_devc.csv', 'Ranz_Marshall_Table_4_10_devc.csv',
    'Ranz_Marshall_Table_4_11_devc.csv', 'Ranz_Marshall_Table_4_12_devc.csv',
    'Ranz_Marshall_Table_4_13_devc.csv'
]
process_table(filenames4, 'Ranz_Marshall_Table_4.csv')
