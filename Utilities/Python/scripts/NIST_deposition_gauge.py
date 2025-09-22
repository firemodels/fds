#!/usr/bin/env python3

# Converted from
# Overholt
# 10-30-2012
# NIST_deposition_gauge.m

import os
import pandas as pd

# Equivalent to Matlab script
def main():
    outdir = os.path.normpath('../../../out/NIST_Deposition_Gauge/') + os.sep

    filenames = [
        'NIST_SDG-2p5-100-p39a_devc.csv', 'NIST_SDG-2p5-100-p39b_devc.csv', 'NIST_SDG-2p5-100-p39c_devc.csv', 'NIST_SDG-2p5-100-p39d_devc.csv',
        'NIST_SDG-5p0-100-p39a_devc.csv', 'NIST_SDG-5p0-100-p39b_devc.csv', 'NIST_SDG-5p0-100-p39c_devc.csv', 'NIST_SDG-5p0-100-p39d_devc.csv',
        'NIST_SDG-10p0-100-p39a_devc.csv', 'NIST_SDG-10p0-100-p39b_devc.csv', 'NIST_SDG-10p0-100-p39c_devc.csv',
        'NIST_SDG-2p5-200-p39a_devc.csv', 'NIST_SDG-2p5-200-p39b_devc.csv', 'NIST_SDG-2p5-200-p39c_devc.csv', 'NIST_SDG-2p5-200-p39d_devc.csv',
        'NIST_SDG-5p0-200-p39a_devc.csv', 'NIST_SDG-5p0-200-p39b_devc.csv', 'NIST_SDG-5p0-200-p39c_devc.csv', 'NIST_SDG-5p0-200-p39d_devc.csv',
        'NIST_SDG-10p0-200-p39a_devc.csv', 'NIST_SDG-10p0-200-p39b_devc.csv', 'NIST_SDG-10p0-200-p39c_devc.csv'
    ]

    # Check for missing files
    skip_case = False
    for fn in filenames:
        fullpath = os.path.join(outdir, fn)
        if not os.path.isfile(fullpath):
            print(f"Error: File {fullpath} does not exist. Skipping case.")
            skip_case = True
    if skip_case:
        return

    # SLPM values
    slpm = [2.5,2.5,2.5,2.5,5,5,5,5,10,10,10,2.5,2.5,2.5,2.5,5,5,5,5,10,10,10]

    # Preallocate “cell arrays”
    sensor1 = [None]*len(filenames)
    sensor2 = [None]*len(filenames)
    sensor3 = [None]*len(filenames)
    sensor4 = [None]*len(filenames)

    # Equivalent of for i=1:22
    for i in range(len(filenames)):
        fullpath = os.path.join(outdir, filenames[i])
        df = pd.read_csv(fullpath, skiprows=2)
        last_row = df.iloc[-1]
        sensor1[i] = last_row.iloc[2]  # Matlab data(end,3)
        sensor2[i] = last_row.iloc[3]
        sensor3[i] = last_row.iloc[4]
        sensor4[i] = last_row.iloc[5]

    # Write out results file (100)
    out100 = os.path.join(outdir, 'SDG_All_Tests_100.csv')
    with open(out100, 'w') as f:
        f.write('SLPM,Sensor1,Sensor2,Sensor3,Sensor4\n')
        for i in range(11):  # Matlab for i=1:11
            f.write(f'{slpm[i]:5.3e},{sensor1[i]:5.3e},{sensor2[i]:5.3e},{sensor3[i]:5.3e},{sensor4[i]:5.3e}\n')

    # Write out results file (200)
    out200 = os.path.join(outdir, 'SDG_All_Tests_200.csv')
    with open(out200, 'w') as f:
        f.write('SLPM,Sensor1,Sensor2,Sensor3,Sensor4\n')
        for i in range(11, 22):  # Matlab for i=1:11 for the second block
            f.write(f'{slpm[i]:5.3e},{sensor1[i]:5.3e},{sensor2[i]:5.3e},{sensor3[i]:5.3e},{sensor4[i]:5.3e}\n')

if __name__ == '__main__':
    main()
