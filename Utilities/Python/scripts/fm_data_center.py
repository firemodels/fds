#!/usr/bin/env python
"""
fm_data_center.py
Jason Floyd
"""

# Generate FDS summary data for scatterplots

import pandas as pd
import numpy as np
import fdsplotlib
import os

outdir = '../../../out/FM_FPRF_Datacenter/'
expdir = '../../../exp/FM_FPRF_Datacenter/'

# Experimental data
exp_data = pd.read_csv(os.path.join(expdir, 'fm_exp.csv'), skiprows=1, header=None).values.flatten()

# --- Low flow test ---
fds_data = pd.read_csv(os.path.join(outdir, 'FM_Datacenter_Veltest_Low_devc.csv'),
                       skiprows=13, header=None).values
n_fds_data = fds_data.shape[0]

# compute average pressures
fds_avg = []
for i in range(1, 7):  # MATLAB 2:7 → Python 1:7
    fds_avg.append(np.mean(fds_data[:, i]))

fds_header = ['Time', 'Low SF-CA', 'Low HA-CP', 'High SF-CA', 'High HA-CP']

fds_out = np.zeros(23)  # preallocate (we fill later)
fds_out[0] = 100
fds_out[1] = fds_avg[0] - fds_avg[1]
fds_out[2] = fds_avg[4] - fds_avg[5]

# --- High flow test ---
fds_data = pd.read_csv(os.path.join(outdir, 'FM_Datacenter_Veltest_High_devc.csv'),
                       skiprows=13, header=None).values
for i in range(1, 7):
    fds_avg[i - 1] = np.mean(fds_data[:, i])

fds_out[3] = fds_avg[0] - fds_avg[1]
fds_out[4] = fds_avg[4] - fds_avg[5]

# Write the first 5 outputs
with open(os.path.join(outdir, 'FM_Datacenter_fds_data.csv'), 'w') as fid:
    fid.write(', '.join(fds_header) + '\n')
    fid.write('{:.0f}, {:.3f}, {:.3f}, {:.3f}, {:.3f}\n'.format(*fds_out[0:5]))

# --- get soot values ---
def mean_soot(filename, skiprows, cols):
    data = pd.read_csv(os.path.join(outdir, filename), skiprows=skiprows, header=None).values
    return [np.mean(data[:, c - 1]) * 1_000_000 for c in cols]  # MATLAB is 1-based

cols = [28, 46, 65]
fds_out[5:8] = mean_soot('FM_Datacenter_Low_C3H6_SF_devc.csv', 15, cols)
fds_out[8:11] = mean_soot('FM_Datacenter_High_C3H6_SF_devc.csv', 15, cols)
fds_out[11:14] = mean_soot('FM_Datacenter_Low_C3H6_HA_devc.csv', 15, cols)
fds_out[14:17] = mean_soot('FM_Datacenter_High_C3H6_HA_devc.csv', 15, cols)
fds_out[17:20] = mean_soot('FM_Datacenter_Low_Cable_SF_devc.csv', 13, cols)
fds_out[20:23] = mean_soot('FM_Datacenter_High_Cable_SF_devc.csv', 13, cols)

# --- plotting ---
x = np.array([0.01, 0.122, 0.2, 0.3, 0.5, 1, 5, 10, 50, 100, 500, 1000])
logx = np.minimum(np.maximum(-2.7, np.log(x + 0.00001)), 1)
errx = 0.01 * (3.8184 * logx ** 2 - 7.7783 * logx + 14.346)
toterr = np.sqrt(errx ** 2 + 0.1 ** 2 + 0.1 ** 2 + 0.05 ** 2)
xerrp = x + 2 * toterr * x
xerrm = np.maximum(0.00001, x - 2 * toterr * x)

fig = fdsplotlib.plot_to_fig(x,x,
                             plot_type='loglog',
                             marker_style='k-',
                             x_min=0.01, x_max=300, y_min=0.01, y_max=300,
                             x_label='Measured Soot Concentration (mg/m$^3$)',
                             y_label='Predicted Soot Concentration (mg/m$^3$)',
                             legend_location='lower right')

fdsplotlib.plot_to_fig(x, xerrp, figure_handle=fig, marker_style='k--')
fdsplotlib.plot_to_fig(x, xerrm, figure_handle=fig, marker_style='k--')
fdsplotlib.plot_to_fig(exp_data[5:8], fds_out[5:8], figure_handle=fig, marker_style='ro', data_label='C3H6 Low SF')
fdsplotlib.plot_to_fig(exp_data[8:11], fds_out[8:11], figure_handle=fig, marker_style='r+', data_label='C3H6 High SF')
fdsplotlib.plot_to_fig(exp_data[11:14], fds_out[11:14], figure_handle=fig, marker_style='bo', data_label='C3H6 Low HA')
fdsplotlib.plot_to_fig(exp_data[14:17], fds_out[14:17], figure_handle=fig, marker_style='b+', data_label='C3H6 High HA')
fdsplotlib.plot_to_fig(exp_data[17:20], fds_out[17:20], figure_handle=fig, marker_style='go', data_label='Cable Low SF')
fdsplotlib.plot_to_fig(exp_data[20:23], fds_out[20:23], figure_handle=fig, marker_style='g+', data_label='Cable High SF')

fig.savefig('../../Manuals/FDS_Validation_Guide/SCRIPT_FIGURES/FM_FPRF_Datacenter/FM_Datacenter_Soot.pdf')

# --- write pressure.tex ---
filename = '../../Manuals/FDS_Validation_Guide/SCRIPT_FIGURES/FM_FPRF_Datacenter/pressure.tex'
with open(filename, 'w') as fid:
    pres_dump = np.zeros(12)
    pres_dump[0] = fds_out[1]
    pres_dump[1] = exp_data[1]
    pres_dump[2] = exp_data[1] * .19
    pres_dump[3] = fds_out[2]
    pres_dump[4] = exp_data[2]
    pres_dump[5] = exp_data[2] * .19
    pres_dump[6] = fds_out[3]
    pres_dump[7] = exp_data[3]
    pres_dump[8] = exp_data[3] * .1
    pres_dump[9] = fds_out[4]
    pres_dump[10] = exp_data[4]
    pres_dump[11] = exp_data[4] * .1

    fid.write('\\begin{center}\n')
    fid.write('\\begin{tabular}{|c|c|c|c|c|} \\hline\n')
    fid.write('Fan Speed & FDS SF to CA & Exp SF to CA & FDS HA to CP & Exp HA to CP \\\\\n')
    fid.write('  & (Pa) & (Pa) & (Pa) & (Pa)  \\\\ \\hline\\hline\n')
    fid.write('78 ACH & {:5.1f} & {:5.1f} $\\pm$ {:5.1f} & {:5.1f} & {:5.1f} $\\pm$ {:5.1f} \\\\\n'.format(*pres_dump[0:6]))
    fid.write('265 ACH & {:5.1f} & {:5.1f} $\\pm$ {:5.1f} & {:5.1f} & {:5.1f} $\\pm$ {:5.1f} \\\\\n'.format(*pres_dump[6:12]))
    fid.write('\\hline\n')
    fid.write('\\end{tabular}\n')
    fid.write('\\end{center}\n')

