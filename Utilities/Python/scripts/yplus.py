
# Create plot using the output of Verification/Turbulence/yplus_N.fds

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import fdsplotlib

plot_style = fdsplotlib.get_plot_style('fds')

outdir = '../../Verification/Turbulence/'
pltdir = '../../Manuals/FDS_Verification_Guide/SCRIPT_FIGURES/'

git_file = outdir + 'yplus_8_git.txt'
version_string = fdsplotlib.get_version_string(git_file)

M_8  = pd.read_csv(os.path.join(outdir, 'yplus_8_devc.csv'),  skiprows=2)
M_16 = pd.read_csv(os.path.join(outdir, 'yplus_16_devc.csv'), skiprows=2)
M_32 = pd.read_csv(os.path.join(outdir, 'yplus_32_devc.csv'), skiprows=2)

yp = np.array([M_8.iloc[-1, 1], M_16.iloc[-1, 1], M_32.iloc[-1, 1]])

n = np.array([8, 16, 32])
y = 0.5 * (1.0 / n)
mu = 1
rho = 1.199
tau_w = 0.5
u_tau = np.sqrt(tau_w / rho)
d_nu = (mu / rho) / u_tau
yp_exact = y / d_nu

fig = fdsplotlib.plot_to_fig(x_data=yp_exact, y_data=yp_exact, marker_style='ko', data_label='Exact',
                             x_min=0.01, x_max=0.05, y_min=0.01, y_max=0.05,
                             revision_label=version_string,
                             x_label='$y^+$ specified',
                             y_label='$y^+$ predicted')

fdsplotlib.plot_to_fig(x_data=yp_exact, y_data=yp, figure_handle=fig, marker_style='k-', data_label='FDS')

plt.savefig(os.path.join(pltdir, 'yplus.pdf'), format='pdf')
plt.close()

