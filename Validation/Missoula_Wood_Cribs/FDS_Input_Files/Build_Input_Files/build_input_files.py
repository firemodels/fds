#6/2-21/2022 Noelle Crump, to convert data from Deep Fuel Bed Exp csv file
#into a paramfile to run through swaps.py to create input files for DFB
#code based on The MATLAB script of the same purpose by Ruddy Mell

import math
import pandas as pd
import os

# Build input files: run swaps.py
os.system('python ../../../../Utilities/Input_File_Tools/swaps.py')

# Move inpupt files up one level to FDS_Input_Files
os.system('mv ./missoula_cribs*.fds ../.')

