
# Make Smokeview images for the FDS Manuals.
# The script gets the smokeview executable appropriate for your OS from the smv repository 
# unless you provide an optional full path to the executable as an argument.

import pandas as pd
import os
import subprocess
import shutil
import platform
import argparse

parser = argparse.ArgumentParser(description="A script to generate Smokeview images")
parser.add_argument("message", nargs="?", default="null", help="Optional smokeview path")
args = parser.parse_args()
smokeview_path = args.message

os_name = platform.system()

if os_name == "Linux":
    if shutil.which("xvfb-run") is None:
        raise FileNotFoundError("xvfb-run is not installed. Please install xvfb package.")

bindir = '../../../smv/Build/for_bundle'
smvdir = '../../../smv/Build/smokeview/'
outdir = '../../Verification/'
original_dir = os.getcwd()

df = pd.read_csv(outdir + 'scripts/FDS_Pictures.csv', header=None)

folder = df[0].values
case = df[1].values

if smokeview_path != "null":
    print("Using "+smokeview_path)
elif os_name == "Linux":
    smokeview_path = smvdir + 'intel_linux_64/smokeview_linux_64'
elif os_name == "Darwin":
    smokeview_path = smvdir + 'gnu_osx_64/smokeview_osx_64'
elif os_name == "Windows":
    smokeview_path = smvdir + 'intel_win_64/smokeview_win_64'

for i in range(len(folder)):
    print('generating smokeview image ' + case[i])
    os.chdir(outdir + folder[i])
    if os_name == "Linux":
        try:
            subprocess.run(['xvfb-run','-w','10','-s','-fp /usr/share/X11/fonts/misc -screen 0 1280x1024x24','-a',smokeview_path,
                       '-bindir',bindir,'-runscript', case[i] ], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        except subprocess.CalledProcessError as e:
            print(f"Smokeview failed with return code {e.returncode}")
        except FileNotFoundError:
            print(f"Smokeview executable not found: {smokeview_path}")
    else:
        subprocess.run([smokeview_path,'-bindir',bindir,'-runscript',case[i]], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    os.chdir(original_dir)

