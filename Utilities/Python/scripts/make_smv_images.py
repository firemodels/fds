
# Make Smokeview images for the FDS Manuals.
# The script uses the command "smokeview" unless you provide an optional name.

import pandas as pd
import os
import subprocess
import shutil
import platform
import argparse

parser = argparse.ArgumentParser(description="A script to generate Smokeview images")
parser.add_argument("message", nargs="?", default="smokeview", help="Optional smokeview name")
args = parser.parse_args()
shell_command = args.message

os_name = platform.system()

if os_name == "Linux":
    if shutil.which("xvfb-run") is None:
        raise FileNotFoundError("xvfb-run is not installed. Please install xvfb package.")

bindir = '../../../smv/Build/for_bundle'
outdir = '../../Verification/'
original_dir = os.getcwd()

df = pd.read_csv(outdir + 'scripts/FDS_Pictures.csv', header=None)

folder = df[0].values
case = df[1].values

if os_name == "Linux":
    result = subprocess.run(["bash", "-i", "-c", f"alias {shell_command} 2>/dev/null | sed -E \"s/alias {shell_command}='(.*)'/\\1/\""],
                            capture_output=True, text=True)
    smokeview_path=result.stdout.strip()
else:
    smokeview_path = shutil.which(shell_command)

for i in range(len(folder)):
    print('generating smokeview image ' + case[i])
    os.chdir(outdir + folder[i])
    if os_name == "Linux":
        try:
            subprocess.run(['xvfb-run','-w','10','-s','-fp /usr/share/X11/fonts/misc -screen 0 1280x1024x24','-a',smokeview_path,
                       '-bindir',bindir,'-runscript', case[i] ], check=True)
        except subprocess.CalledProcessError as e:
            print(f"Smokeview failed with return code {e.returncode}")
        except FileNotFoundError:
            print(f"Smokeview executable not found: {smokeview_path}")

    else:
        subprocess.run([smokeview_path,'-bindir',bindir,'-runscript',case[i]], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    os.chdir(original_dir)

