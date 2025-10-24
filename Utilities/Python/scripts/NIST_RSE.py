
# Process FDS simulation data of NIST RSE 1994 Experiments

import numpy as np
import pandas as pd
import os

outdir = '../../../out/NIST_RSE_1994/'

HRR = [50, 75, 100, 150, 200, 300, 400, 500, 600]

# Initialize 3D arrays for storing data
# Will determine dimensions from first file read
num_cols = None

#---50kW---%
RSE_50_struct = pd.read_csv(outdir + 'NIST_RSE_1994_50_RI=5_devc.csv') # FDS 50 kW
num_cols = RSE_50_struct.shape[1]
RSE_data_RI_5 = np.zeros((1, num_cols, 9))
RSE_data_RI_10 = np.zeros((1, num_cols, 9))
RSE_data_RI_20 = np.zeros((1, num_cols, 9))
RSE_data_RI_5[:, :, 0] = RSE_50_struct.iloc[-1, :].values.reshape(1, -1) # data

RSE_50_struct = pd.read_csv(outdir + 'NIST_RSE_1994_50_RI=10_devc.csv') # FDS 50 kW
RSE_data_RI_10[:, :, 0] = RSE_50_struct.iloc[-1, :].values.reshape(1, -1) # data

RSE_50_struct = pd.read_csv(outdir + 'NIST_RSE_1994_50_RI=20_devc.csv') # FDS 50 kW
RSE_data_RI_20[:, :, 0] = RSE_50_struct.iloc[-1, :].values.reshape(1, -1) # data

#---75kW---%
RSE_75_struct = pd.read_csv(outdir + 'NIST_RSE_1994_75_RI=5_devc.csv') # FDS 75 kW
RSE_data_RI_5[:, :, 1] = RSE_75_struct.iloc[-1, :].values.reshape(1, -1) # data

RSE_75_struct = pd.read_csv(outdir + 'NIST_RSE_1994_75_RI=10_devc.csv') # FDS 75 kW
RSE_data_RI_10[:, :, 1] = RSE_75_struct.iloc[-1, :].values.reshape(1, -1) # data

RSE_75_struct = pd.read_csv(outdir + 'NIST_RSE_1994_75_RI=20_devc.csv') # FDS 75 kW
RSE_data_RI_20[:, :, 1] = RSE_75_struct.iloc[-1, :].values.reshape(1, -1) # data

#---100kW---%
RSE_100_struct = pd.read_csv(outdir + 'NIST_RSE_1994_100_RI=5_devc.csv') # FDS 100 kW
RSE_data_RI_5[:, :, 2] = RSE_100_struct.iloc[-1, :].values.reshape(1, -1) # data

RSE_100_struct = pd.read_csv(outdir + 'NIST_RSE_1994_100_RI=10_devc.csv') # FDS 100 kW
RSE_data_RI_10[:, :, 2] = RSE_100_struct.iloc[-1, :].values.reshape(1, -1) # data

RSE_100_struct = pd.read_csv(outdir + 'NIST_RSE_1994_100_RI=20_devc.csv') # FDS 100 kW
RSE_data_RI_20[:, :, 2] = RSE_100_struct.iloc[-1, :].values.reshape(1, -1) # data

#---150kW---%
RSE_150_struct = pd.read_csv(outdir + 'NIST_RSE_1994_150_RI=5_devc.csv') # FDS 150 kW
RSE_data_RI_5[:, :, 3] = RSE_150_struct.iloc[-1, :].values.reshape(1, -1) # data

RSE_150_struct = pd.read_csv(outdir + 'NIST_RSE_1994_150_RI=10_devc.csv') # FDS 150 kW
RSE_data_RI_10[:, :, 3] = RSE_150_struct.iloc[-1, :].values.reshape(1, -1) # data

RSE_150_struct = pd.read_csv(outdir + 'NIST_RSE_1994_150_RI=20_devc.csv') # FDS 150 kW
RSE_data_RI_20[:, :, 3] = RSE_150_struct.iloc[-1, :].values.reshape(1, -1) # data

#---200kW---%
RSE_200_struct = pd.read_csv(outdir + 'NIST_RSE_1994_200_RI=5_devc.csv') # FDS 200 kW
RSE_data_RI_5[:, :, 4] = RSE_200_struct.iloc[-1, :].values.reshape(1, -1) # data

RSE_200_struct = pd.read_csv(outdir + 'NIST_RSE_1994_200_RI=10_devc.csv') # FDS 200 kW
RSE_data_RI_10[:, :, 4] = RSE_200_struct.iloc[-1, :].values.reshape(1, -1) # data

RSE_200_struct = pd.read_csv(outdir + 'NIST_RSE_1994_200_RI=20_devc.csv') # FDS 200 kW
RSE_data_RI_20[:, :, 4] = RSE_200_struct.iloc[-1, :].values.reshape(1, -1) # data

#---300kW---%
RSE_300_struct = pd.read_csv(outdir + 'NIST_RSE_1994_300_RI=5_devc.csv') # FDS 300 kW
RSE_data_RI_5[:, :, 5] = RSE_300_struct.iloc[-1, :].values.reshape(1, -1) # data

RSE_300_struct = pd.read_csv(outdir + 'NIST_RSE_1994_300_RI=10_devc.csv') # FDS 300 kW
RSE_data_RI_10[:, :, 5] = RSE_300_struct.iloc[-1, :].values.reshape(1, -1) # data

RSE_300_struct = pd.read_csv(outdir + 'NIST_RSE_1994_300_RI=20_devc.csv') # FDS 300 kW
RSE_data_RI_20[:, :, 5] = RSE_300_struct.iloc[-1, :].values.reshape(1, -1) # data

#---400kW---%
RSE_400_struct = pd.read_csv(outdir + 'NIST_RSE_1994_400_RI=5_devc.csv') # FDS 400 kW
RSE_data_RI_5[:, :, 6] = RSE_400_struct.iloc[-1, :].values.reshape(1, -1) # data

RSE_400_struct = pd.read_csv(outdir + 'NIST_RSE_1994_400_RI=10_devc.csv') # FDS 400 kW
RSE_data_RI_10[:, :, 6] = RSE_400_struct.iloc[-1, :].values.reshape(1, -1) # data

RSE_400_struct = pd.read_csv(outdir + 'NIST_RSE_1994_400_RI=20_devc.csv') # FDS 400 kW
RSE_data_RI_20[:, :, 6] = RSE_400_struct.iloc[-1, :].values.reshape(1, -1) # data

#---500kW---%
RSE_500_struct = pd.read_csv(outdir + 'NIST_RSE_1994_500_RI=5_devc.csv') # FDS 500 kW
RSE_data_RI_5[:, :, 7] = RSE_500_struct.iloc[-1, :].values.reshape(1, -1) # data

RSE_500_struct = pd.read_csv(outdir + 'NIST_RSE_1994_500_RI=10_devc.csv') # FDS 500 kW
RSE_data_RI_10[:, :, 7] = RSE_500_struct.iloc[-1, :].values.reshape(1, -1) # data

RSE_500_struct = pd.read_csv(outdir + 'NIST_RSE_1994_500_RI=20_devc.csv') # FDS 500 kW
RSE_data_RI_20[:, :, 7] = RSE_500_struct.iloc[-1, :].values.reshape(1, -1) # data

#---600kW---%
RSE_600_struct = pd.read_csv(outdir + 'NIST_RSE_1994_600_RI=5_devc.csv') # FDS 600 kW
RSE_data_RI_5[:, :, 8] = RSE_600_struct.iloc[-1, :].values.reshape(1, -1) # data

RSE_600_struct = pd.read_csv(outdir + 'NIST_RSE_1994_600_RI=10_devc.csv') # FDS 600 kW
RSE_data_RI_10[:, :, 8] = RSE_600_struct.iloc[-1, :].values.reshape(1, -1) # data

RSE_600_struct = pd.read_csv(outdir + 'NIST_RSE_1994_600_RI=20_devc.csv') # FDS 600 kW
RSE_data_RI_20[:, :, 8] = RSE_600_struct.iloc[-1, :].values.reshape(1, -1) # data

#------------------------------------------------
# Building FDS data file with average species/temp
# Function of HRR
#------------------------------------------------

# Calculate number of rows needed for RSE_Results
# Row 1: HRR values
# Rows 2-16:  RI_5  data (15 rows, columns 2-16 from original data)
# Rows 17-31: RI_10 data (15 rows, columns 2-16 from original data)
# Rows 32-46: RI_20 data (15 rows, columns 2-16 from original data)
num_data_rows = 2 * (num_cols - 1) + num_cols
RSE_Results = np.zeros((num_data_rows, 9))

for i in range(9):
    RSE_Results[0, i] = HRR[i]
    for j in range(1, 16): 
        RSE_Results[j, i] = RSE_data_RI_5[0, j, i]
        RSE_Results[num_cols + j - 1, i] = RSE_data_RI_10[0, j, i]
        RSE_Results[2 * (num_cols - 1) + j, i] = RSE_data_RI_20[0, j, i]

header1 = ['HRR', 'O2Rear_FDS_RI_5', 'CO2Rear_FDS_RI_5', 'CORear_FDS_RI_5', 'UHRear_FDS_RI_5', 'H2ORear_FDS_RI_5',
    'O2Front_FDS_RI_5', 'CO2Front_FDS_RI_5', 'COFront_FDS_RI_5', 'UHFront_FDS_RI_5', 'H2OFront_FDS_RI_5',
    'TRSampA_FDS_RI_5', 'TRSampBB_FDS_RI_5', 'TFSampA_FDS_RI_5', 'TFSampBB_FDS_RI_5', 'ITER_RI=5', 'O2Rear_FDS_RI_10',
    'CO2Rear_FDS_RI_10', 'CORear_FDS_RI_10', 'UHRear_FDS_RI_10', 'H2ORear_FDS_RI_10', 'O2Front_FDS_RI_10',
    'CO2Front_FDS_RI_10', 'COFront_FDS_RI_10', 'UHFront_FDS_RI_10', 'H2OFront_FDS_RI_10', 'TRSampA_FDS_RI_10',
    'TRSampBB_FDS_RI_10', 'TFSampA_FDS_RI_10', 'TFSampBB_FDS_RI_10', 'ITER_RI=10', 'O2Rear_FDS_RI_20', 'CO2Rear_FDS_RI_20',
    'CORear_FDS_RI_20', 'UHRear_FDS_RI_20', 'H2ORear_FDS_RI_20', 'O2Front_FDS_RI_20', 'CO2Front_FDS_RI_20',
    'COFront_FDS_RI_20', 'UHFront_FDS_RI_20', 'H2OFront_FDS_RI_20', 'TRSampA_FDS_RI_20', 'TRSampBB_FDS_RI_20',
    'TFSampA_FDS_RI_20', 'TFSampBB_FDS_RI_20', 'ITER_RI=20']

filename1 = outdir + 'NIST_RSE_1994_FDS.csv'
fid = open(filename1, 'w')
fid.write(', '.join(header1) + '\n')
for i in range(9):
    row_values = ['{:f}'.format(val) for val in RSE_Results[:, i]]
    fid.write(', '.join(row_values) + '\n ')
fid.close()

