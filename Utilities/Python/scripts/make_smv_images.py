
# Make Smokeview images for the FDS Manuals

import pandas as pd
import os
import subprocess
import shutil
import platform

os_name = platform.system()

if os_name == "Linux":
    if shutil.which("xvfb-run") is None:
        raise FileNotFoundError("xvfb-run is not installed. Please install xvfb package.")

outdir = '../../Verification/'
original_dir = os.getcwd()

df = pd.read_csv(outdir + 'scripts/FDS_Pictures.csv', header=None)

folder = df[0].values
case = df[1].values

smokeview_path = shutil.which('smokeview')

for i in range(len(folder)):
    os.chdir(outdir + folder[i])
    if os_name == "Linux":
        subprocess.run(['xvfb-run',smokeview_path,'-runscript',case[i]])
    else:
        subprocess.run([smokeview_path,'-runscript',case[i]])
    os.chdir(original_dir)

