#!/usr/bin/env python3
"""
swaps.py (SWAP Strings)
by Randy McDermott
June 2021

Perform string replacements in an FDS input file (fdsbasefile) based on a parameter file (paramfile).

CAUTION: The string replacements must have *unique* parameter names.  For example, do not simply use
         an ordered set of parameters with, say, param1 going up to param10, because swaps thinks that the
         first 6 elements of param10 are param1 and gives the wrong result.
"""

import pandas as pd
import shutil

# change the name of "paramfile.csv"
paramfile = 'paramfile.csv'

df = pd.read_csv(paramfile, sep=' *, *', engine='python', dtype=str)

fdsbasefile = df.columns[0]
# print(fdsbasefile)

for irow in df.index:

    fdsoutfile = df.loc[irow,fdsbasefile]
    print('Generating '+fdsoutfile)
    # copy the base file to a new file
    shutil.copy(fdsbasefile,fdsoutfile)

    for param in df.columns[1:]:
        # Read the content of file
        f1 = open(fdsoutfile, 'r')
        input_data = f1.read()
        f1.close()
        # Replace the target string
        input_data = input_data.replace(param, df.loc[irow,param])
        # Write the output to same file
        f2 = open(fdsoutfile, 'w')
        f2.write(input_data)
        f2.close()
