# N. Crump 6/11/2021 - 8/5/2021
# Askervein Hill case file processor (Python version)
#
# Script Process:
#-Read in the decv.csv file from the out directory
#-Build the output matrix with sorted rows
#-Insert the Distance from Hiltop column
#-output new .csv file

# Navigation
#-10m Post lines A, AA, B
     #ASW-ANE:   SPEED DIRECTION UPWASH SIGM_SPEED SIGMw
     #BNW-BSE:   SPEED
     #AASW-AANE: SPEED
#-Mean flow runs MF and TU
     #TU:   SPEED DIRECTION UPWASH SIGM_SPEED SIGMw
     #MF:   SPEED
#-50m Profile at HT 
   #     Data:   SPEED, SIGM_SPEED
#-FRG 17m Profile at CP' 
     #TU:   SPEED DIRECTION UPWASH SIGM_SPEED
     #MF:   SPEED
#-UK 30m Profile at ASW60 
   #     Data:   SPEED DIRECTION UPWASH SIGM_SPEED

# Questions Checklist:
# reading in the column headers is very strange '"..."' the code will
# not work without the extra "". To see this, >> A.colheaders
 
# Debugging notes:
# -The most common error will be dimensions of arrays cannot be concacted
# this is most often caused by a failure to read in the the script. 
#   Check if Time_Stop is set to the correct value
#   Check if the file has been saved in column veiwing format
#   Check if the colheaders have '".."' or '..' with >> A.colheaders
#   Check values with >> whos DIST DIRECTION UPWASH AVG_SPEED SIGM_SPEED SIGMw

import os
import numpy as np
import pandas as pd

outdir = '../../../out/Askervein_Hill/'

# --- read CSV with two header lines like MATLAB ---
filename = 'Askervein_TU03A_16m_devc.csv'
fpath = os.path.join(outdir, filename)

# First grab header line with column names
with open(fpath, 'r') as f:
    header_lines = [next(f) for _ in range(2)]
colheaders = header_lines[1].strip().split(',')

# Then read the data only
A = pd.read_csv(fpath, skiprows=2, header=None)

# replicate MATLAB variable ID
ID = filename.replace('Askervein_TU03A_', '').replace('_devc.csv', '')

# time slot selection
Time_idx = colheaders.index('Time')
Times = A.iloc[:, Time_idx].values
TimeStop = Times.max()
slot = np.where((Times > TimeStop - 1) & (Times < TimeStop + 1))[0]

# ========== A-Line 10m ==========
NewName = f'ASW-ANE_TU03A_{ID}.csv'
DIST = np.array(["(m)", "DIST", -85, -60, -50, -35, -20, -10, 0, 10, 20, 40], dtype=object)
TOWERS = np.array(["ASW85","ASW UK 30 m tower 10m","ASW50","ASW35","ASW20","ASW10","HT 10 m t","ANE10","ANE20","ANE40"], dtype=object)
QUANTITY = np.array(["SPEED ","SIGM SPEED ","SIGMw ","AVGu ","AVGv ","AVGw "], dtype=object)

# build “colls” array like MATLAB
colls = []
f = '"'
for this_quant in QUANTITY:
    new_row = []
    for this_tower in TOWERS:
        this_name = f + this_quant + this_tower + f
        idx = colheaders.index(this_name)
        value = A.iloc[slot, idx].values[0]  # single time slot
        new_row.append(value)
    colls.append(new_row)
colls = np.array(colls)

AVG_SPEED = np.concatenate((["(m/s)","SPEED"], colls[0,:]))
SIGM_SPEED = np.concatenate((["(m/s)","SIGM_SPEED"], colls[1,:]))
SIGMw = np.concatenate((["(m/s)","SIGMw"], colls[2,:]))
AVGu = colls[3,:].astype(float)
AVGv = colls[4,:].astype(float)
AVGw = colls[5,:].astype(float)

# secondary variables
HORIZ_VELO = np.sqrt(AVGv**2 + AVGu**2)
UPWASH = np.concatenate((["(deg)","UPWASH"], np.degrees(np.arctan(AVGw/HORIZ_VELO))))
DIRECTION_S0 = np.sign(AVGu)*90 + 180 - np.degrees(np.arctan(AVGv/AVGu))
DIRECTION = np.concatenate((["(deg)","DIRECTION"], DIRECTION_S0))

# combine columns
FINAL = pd.DataFrame({
    'DIST': DIST,
    'DIRECTION': DIRECTION,
    'UPWASH': UPWASH,
    'SPEED': AVG_SPEED,
    'SIGM_SPEED': SIGM_SPEED,
    'SIGMw': SIGMw
})

FINAL.to_csv(os.path.join(outdir, NewName), index=False, header=False)
