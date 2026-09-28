
# Make Smokeview images for the FDS Manuals.
# The script gets the smokeview executable appropriate for your OS from the smv repository 
# unless you provide an optional full path to the executable as an argument.
# After images are created, a comparison is made with the reference figures in the fig repository.

import pandas as pd
import os
import subprocess
import shutil
import platform
import argparse
import sys
from PIL import Image
import numpy as np

def compare_images(image1_path, image2_path, threshold=10):
    img1 = Image.open(image1_path).convert("RGB")
    img2 = Image.open(image2_path).convert("RGB")

    if img1.size != img2.size:
        raise ValueError("Images must have the same dimensions")

    arr1 = np.array(img1, dtype=np.int16)
    arr2 = np.array(img2, dtype=np.int16)

    diff = np.abs(arr1 - arr2)
    pixel_diff = np.max(diff, axis=2)

    changed_pixels = np.sum(pixel_diff > threshold)
    total_pixels = pixel_diff.size

    difference_percentage = (changed_pixels / total_pixels) * 100

    mse = np.mean((arr1 - arr2) ** 2)

    return {
        "mse": mse,
        "difference_percentage": difference_percentage
    }


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
    smokeview_path = smvdir + 'intel_linux/smokeview_linux'
elif os_name == "Darwin":
    smokeview_path = smvdir + 'gnu_osx/smokeview_osx'
elif os_name == "Windows":
    smokeview_path = smvdir + 'intel_win/smokeview_win'

for i in range(len(folder)):
    print('generating smokeview image ' + case[i])
    os.chdir(outdir + folder[i])
    try:
        if os_name == "Linux":
            subprocess.run(['xvfb-run','-w','2','-s','-fp /usr/share/X11/fonts/misc -screen 0 1280x1024x24','-a',smokeview_path,
                    '-bindir',bindir,'-runscript', '-render_overwrite', case[i] ], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        else:
            subprocess.run([smokeview_path,'-bindir',bindir,'-runscript','-render_overwrite', case[i]], 
                           check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    except subprocess.CalledProcessError as e:
        print(f"Error: Smokeview failed with return code {e.returncode}")
    except FileNotFoundError:
        print(f"Error: Smokeview executable not found: {smokeview_path}")

    os.chdir(original_dir)

# Compare images just created against the reference images in the fig repository

from pathlib import Path

mandir = '../../Manuals/'
refdir = '../../../fig/fds/Reference_Figures/'
directories = [mandir+'FDS_User_Guide/SCRIPT_FIGURES/' , mandir+'FDS_Verification_Guide/SCRIPT_FIGURES/']

print('comparing images with reference figures...')

for directory in directories:
    dir_path = Path(directory)
    
    for png_file in dir_path.glob("*.png"):

        try:
            results = compare_images(directory+png_file.name, refdir+png_file.name)
            if results['difference_percentage'] > 1:
                print('Warning: ',png_file.name,' changed ',f"{results['difference_percentage']:.2f}%")
        except FileNotFoundError as e:
            print(f"Error: Could not find image file - {e}")
            sys.exit(1)
        except Exception as e:
            print(f"Error: {e}")
            sys.exit(1)

