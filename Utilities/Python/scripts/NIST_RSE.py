"""
------------------------------------------------
C Weinschenk
Process FDS simulation data of NIST RSE 1994 Experiments
Converted from MATLAB
------------------------------------------------
"""

import pandas as pd
import numpy as np

# ------------------------------------------------
# Adding paths for data files
# and post processing scripts
# ------------------------------------------------

outdir = '../../../out/NIST_RSE_1994/'

# ------------------------------------------------
# Read in FDS output files
# ------------------------------------------------

HRR = [50, 75, 100, 150, 200, 300, 400, 500, 600]

# Pre-allocate 3-D arrays like MATLAB (list of lists of lists)
RSE_data_RI_5 = []
RSE_data_RI_10 = []
RSE_data_RI_20 = []

# Helper to read one CSV, return last row (like MATLAB .data(end,:))
def last_row_csv(path):
    df = pd.read_csv(path)
    return df.iloc[-1].to_numpy()

# Each HRR corresponds to a new slice (third dimension in MATLAB)
for i, hrr in enumerate(HRR):
    # RI=5
    data5 = last_row_csv(f"{outdir}NIST_RSE_1994_{hrr}_RI=5_devc.csv")
    RSE_data_RI_5.append(data5)
    # RI=10
    data10 = last_row_csv(f"{outdir}NIST_RSE_1994_{hrr}_RI=10_devc.csv")
    RSE_data_RI_10.append(data10)
    # RI=20
    data20 = last_row_csv(f"{outdir}NIST_RSE_1994_{hrr}_RI=20_devc.csv")
    RSE_data_RI_20.append(data20)

# Convert to NumPy arrays with shape (n_HRR, n_columns)
RSE_data_RI_5 = np.array(RSE_data_RI_5)
RSE_data_RI_10 = np.array(RSE_data_RI_10)
RSE_data_RI_20 = np.array(RSE_data_RI_20)

# ------------------------------------------------
# Building FDS data file with average species/temp
# Function of HRR
# ------------------------------------------------

# Build RSE_Results in a MATLAB-like way
# (This part is a straight translation; adjust if column count differs)
RSE_Results = []

for i in range(len(HRR)):
    col = []
    col.append(HRR[i])  # first entry HRR
    # Loop over remaining columns, stacking RI=5/10/20 values
    for j in range(1, 16):  # j=2:16 in MATLAB → range(1,16) in Python
        col.append(RSE_data_RI_5[i][j])
        col.append(RSE_data_RI_10[i][j])
        col.append(RSE_data_RI_20[i][j])
    RSE_Results.append(col)

# Convert to DataFrame
header1 = ['HRR','O2Rear_FDS_RI_5','CO2Rear_FDS_RI_5','CORear_FDS_RI_5','UHRear_FDS_RI_5','H2ORear_FDS_RI_5',
    'O2Front_FDS_RI_5','CO2Front_FDS_RI_5','COFront_FDS_RI_5','UHFront_FDS_RI_5','H2OFront_FDS_RI_5',
    'TRSampA_FDS_RI_5','TRSampBB_FDS_RI_5','TFSampA_FDS_RI_5','TFSampBB_FDS_RI_5','ITER_RI=5','O2Rear_FDS_RI_10',
    'CO2Rear_FDS_RI_10','CORear_FDS_RI_10','UHRear_FDS_RI_10','H2ORear_FDS_RI_10','O2Front_FDS_RI_10',
    'CO2Front_FDS_RI_10','COFront_FDS_RI_10','UHFront_FDS_RI_10','H2OFront_FDS_RI_10','TRSampA_FDS_RI_10',
    'TRSampBB_FDS_RI_10','TFSampA_FDS_RI_10','TFSampBB_FDS_RI_10','ITER_RI=10','O2Rear_FDS_RI_20','CO2Rear_FDS_RI_20',
    'CORear_FDS_RI_20','UHRear_FDS_RI_20','H2ORear_FDS_RI_20','O2Front_FDS_RI_20','CO2Front_FDS_RI_20',
    'COFront_FDS_RI_20','UHFront_FDS_RI_20','H2OFront_FDS_RI_20','TRSampA_FDS_RI_20','TRSampBB_FDS_RI_20',
    'TFSampA_FDS_RI_20','TFSampBB_FDS_RI_20','ITER_RI=20']

df_out = pd.DataFrame(RSE_Results, columns=header1[:len(RSE_Results[0])])

# ------------------------------------------------
# FDS Data CSV File
# ------------------------------------------------

filename1 = f"{outdir}NIST_RSE_1994_FDS.csv"
df_out.to_csv(filename1, index=False)

