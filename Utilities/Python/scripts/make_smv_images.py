
# Make Smokeview images for the FDS Manuals

import pandas as pd
import os
import subprocess

outdir = '../../Verification/'
original_dir = os.getcwd()

df = pd.read_csv(outdir + 'scripts/FDS_Pictures.csv', header=None)

folder = df[0].values
case = df[1].values

for i in range(len(folder)):
    os.chdir(outdir + folder[i])
    subprocess.run(['smokeview','-runscript',case[i]])
    os.chdir(original_dir)

