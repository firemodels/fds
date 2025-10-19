# Floyd
# 5-23-2014
# rms_cov_corr.m
#
# Converted by Floyd
# 10/16/2025

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import fdsplotlib
import os


datadir = '../../Verification/Controls/'
plotdir = '../../Manuals/FDS_Verification_Guide/SCRIPT_FIGURES/'
filename = datadir + 'rms_cov_corr_devc.csv'

skip_case = False

if not os.path.exists(filename):
   skip_case = True
   print('Error: File ', filename, ' does not exist. Skipping case.')

if skip_case: quit()

df = pd.read_csv(filename, skiprows=2, header=None)

start_index = int(500 / 0.02)
end_index = df.index[-1]

data_slice = df.iloc[start_index:end_index + 1]
u = data_slice[1].values
w = data_slice[2].values

umean = np.mean(u)
wmean = np.mean(w)

udiff = u - umean
wdiff = w - wmean

# Calculate RMS, Covariance, and Correlation
urms = np.sqrt(np.mean(udiff**2))
wrms = np.sqrt(np.mean(wdiff**2))

uwcov = np.mean(udiff * wdiff)

uwcorr = uwcov / (urms * wrms)

urms_fds = df.iloc[end_index, 3]
uwcov_fds = df.iloc[end_index, 4]
uwcorr_fds = df.iloc[end_index, 5]

header='Time,urms,uwcov,uwcorr'
outdata=np.zeros((2,4))
outdata[0,0]=0
outdata[1,0]=1000
outdata[:,1] = urms
outdata[:,2] = uwcov
outdata[:,3] = uwcorr

print(header)
print(outdata)

with open(datadir+'rms_cov_corr.csv','w') as fid:
   fid.write(f"{header}\n")
   string_with_separator = np.array2string(outdata[0,:], separator=', ')[1:-1] # [1:-1] removes brackets
   print(string_with_separator)
   fid.write(f"{string_with_separator}\n")
   string_with_separator = np.array2string(outdata[1,:], separator=', ')[1:-1] # [1:-1] removes brackets
   print(string_with_separator)
   fid.write(f"{string_with_separator}\n")


