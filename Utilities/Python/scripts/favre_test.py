
# compare TEMPORAL_STATISTIC='FAVRE AVERAGE' with brute force Favre average

import pandas as pd
import numpy as np

outdir = '../../Verification/Species/'

# gather Favre averages from line file

L = pd.read_csv(outdir + 'favre_test_line.csv', skiprows=1)
M = pd.read_csv(outdir + 'favre_test_devc.csv', skiprows=1)
t_stats_start = M.iloc[-1, 0] / 2
t = M[M.iloc[:, 0] > t_stats_start].iloc[:, 0].values
T = t[-1] - t[0]

# cell data

YL1 = L.iloc[0, L.columns.get_loc('YTILDE_O2')]
YL2 = L.iloc[1, L.columns.get_loc('YTILDE_O2')]
YL3 = L.iloc[2, L.columns.get_loc('YTILDE_O2')]

YRMSL1 = L.iloc[0, L.columns.get_loc('YO2_RMS')]
YRMSL2 = L.iloc[1, L.columns.get_loc('YO2_RMS')]
YRMSL3 = L.iloc[2, L.columns.get_loc('YO2_RMS')]

RHO1 = M[M.iloc[:, 0] > t_stats_start]['RHO_1'].values
RHO2 = M[M.iloc[:, 0] > t_stats_start]['RHO_2'].values
RHO3 = M[M.iloc[:, 0] > t_stats_start]['RHO_3'].values
RHOYO2_1 = M[M.iloc[:, 0] > t_stats_start]['RHOYO2_1'].values
RHOYO2_2 = M[M.iloc[:, 0] > t_stats_start]['RHOYO2_2'].values
RHOYO2_3 = M[M.iloc[:, 0] > t_stats_start]['RHOYO2_3'].values

# brute force integration for means

NUM1 = 0
NUM2 = 0
NUM3 = 0
DENOM1 = 0
DENOM2 = 0
DENOM3 = 0
for i in range(1, len(t)):
    dt = t[i] - t[i-1]
    NUM1 = NUM1 + RHOYO2_1[i] * dt
    NUM2 = NUM2 + RHOYO2_2[i] * dt
    NUM3 = NUM3 + RHOYO2_3[i] * dt
    DENOM1 = DENOM1 + RHO1[i] * dt
    DENOM2 = DENOM2 + RHO2[i] * dt
    DENOM3 = DENOM3 + RHO3[i] * dt

YTILDE1 = NUM1 / DENOM1
YTILDE2 = NUM2 / DENOM2
YTILDE3 = NUM3 / DENOM3

# compute error and report if necessary

e1 = abs(YL1 - YTILDE1)
e2 = abs(YL2 - YTILDE2)
e3 = abs(YL3 - YTILDE3)

tol = 1e-4
if e1 > tol:
    print(f'Matlab Warning: e1 = {e1} in Species/favre_test mean cell data')
if e2 > tol:
    print(f'Matlab Warning: e2 = {e2} in Species/favre_test mean cell data')
if e3 > tol:
    print(f'Matlab Warning: e3 = {e3} in Species/favre_test mean cell data')

# brute force integration for rms

NUM1 = 0
NUM2 = 0
NUM3 = 0
for i in range(1, len(t)):
    dt = t[i] - t[i-1]
    NUM1 = NUM1 + (RHOYO2_1[i] / RHO1[i] - YTILDE1)**2 * dt
    NUM2 = NUM2 + (RHOYO2_2[i] / RHO2[i] - YTILDE2)**2 * dt
    NUM3 = NUM3 + (RHOYO2_3[i] / RHO3[i] - YTILDE3)**2 * dt

YRMS1 = np.sqrt(NUM1 / T)
YRMS2 = np.sqrt(NUM2 / T)
YRMS3 = np.sqrt(NUM3 / T)

# compute error and report if necessary

e1 = abs(YRMSL1 - YRMS1)
e2 = abs(YRMSL2 - YRMS2)
e3 = abs(YRMSL3 - YRMS3)

tol = 1e-2
if e1 > tol:
    print(f'Matlab Warning: e1 = {e1} in Species/favre_test rms cell data')
if e2 > tol:
    print(f'Matlab Warning: e2 = {e2} in Species/favre_test rms cell data')
if e3 > tol:
    print(f'Matlab Warning: e3 = {e3} in Species/favre_test rms cell data')

# interpolated data

YL1 = L.iloc[0, L.columns.get_loc('YTILDE_O2_INT')]
YL2 = L.iloc[1, L.columns.get_loc('YTILDE_O2_INT')]
YL3 = L.iloc[2, L.columns.get_loc('YTILDE_O2_INT')]

YRMSL1 = L.iloc[0, L.columns.get_loc('YO2_RMS_INT')]
YRMSL2 = L.iloc[1, L.columns.get_loc('YO2_RMS_INT')]
YRMSL3 = L.iloc[2, L.columns.get_loc('YO2_RMS_INT')]

RHO1 = M[M.iloc[:, 0] > t_stats_start]['RHO_1_INT'].values
RHO2 = M[M.iloc[:, 0] > t_stats_start]['RHO_2_INT'].values
RHO3 = M[M.iloc[:, 0] > t_stats_start]['RHO_3_INT'].values
RHOYO2_1 = M[M.iloc[:, 0] > t_stats_start]['RHOYO2_1_INT'].values
RHOYO2_2 = M[M.iloc[:, 0] > t_stats_start]['RHOYO2_2_INT'].values
RHOYO2_3 = M[M.iloc[:, 0] > t_stats_start]['RHOYO2_3_INT'].values

# brute force integration

NUM1 = 0
NUM2 = 0
NUM3 = 0
DENOM1 = 0
DENOM2 = 0
DENOM3 = 0
for i in range(1, len(t)):
    dt = t[i] - t[i-1]
    NUM1 = NUM1 + RHOYO2_1[i] * dt
    NUM2 = NUM2 + RHOYO2_2[i] * dt
    NUM3 = NUM3 + RHOYO2_3[i] * dt
    DENOM1 = DENOM1 + RHO1[i] * dt
    DENOM2 = DENOM2 + RHO2[i] * dt
    DENOM3 = DENOM3 + RHO3[i] * dt

YTILDE1 = NUM1 / DENOM1
YTILDE2 = NUM2 / DENOM2
YTILDE3 = NUM3 / DENOM3

# compute error and report if necessary

e1 = abs(YL1 - YTILDE1)
e2 = abs(YL2 - YTILDE2)
e3 = abs(YL3 - YTILDE3)

tol = 1e-4
if e1 > tol:
    print(f'Matlab Warning: e1 = {e1} in Species/favre_test mean interpolated data')
if e2 > tol:
    print(f'Matlab Warning: e2 = {e2} in Species/favre_test mean interpolated data')
if e3 > tol:
    print(f'Matlab Warning: e3 = {e3} in Species/favre_test mean interpolated data')

# brute force integration for rms

NUM1 = 0
NUM2 = 0
NUM3 = 0
for i in range(1, len(t)):
    dt = t[i] - t[i-1]
    NUM1 = NUM1 + (RHOYO2_1[i] / RHO1[i] - YTILDE1)**2 * dt
    NUM2 = NUM2 + (RHOYO2_2[i] / RHO2[i] - YTILDE2)**2 * dt
    NUM3 = NUM3 + (RHOYO2_3[i] / RHO3[i] - YTILDE3)**2 * dt

YRMS1 = np.sqrt(NUM1 / T)
YRMS2 = np.sqrt(NUM2 / T)
YRMS3 = np.sqrt(NUM3 / T)

# compute error and report if necessary

e1 = abs(YRMSL1 - YRMS1)
e2 = abs(YRMSL2 - YRMS2)
e3 = abs(YRMSL3 - YRMS3)

tol = 1e-2
if e1 > tol:
    print(f'Matlab Warning: e1 = {e1} in Species/favre_test rms interpolated data')
if e2 > tol:
    print(f'Matlab Warning: e2 = {e2} in Species/favre_test rms interpolated data')
if e3 > tol:
    print(f'Matlab Warning: e3 = {e3} in Species/favre_test rms interpolated data')

