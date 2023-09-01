#!/usr/bin/env python3
"""
swaps.py (SWAP Strings)
by Randy McDermott
June 2021

Perform string replacements in an FDS input file (fdsbasefile) based on a parameter file (paramfile).
"""

import pandas as pd
import shutil

# change the name of "paramfile.csv"
paramfile = 'paramfile.csv'

df = pd.read_csv(paramfile, sep=' *, *', engine='python', dtype=str)

# assumes Template.fds is the first column name
fdsbasefile = df.columns[0]
# fdsbasefile = 'Frederica.fds' # or just directly assign the filename
# print(fdsbasefile)

# sort dataframe column order by descending column name length to avoid swapping substrings
cols = sorted(df.columns.tolist(),key=len,reverse=True)
df = df[cols]

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
