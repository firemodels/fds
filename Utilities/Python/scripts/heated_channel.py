
# Create plots for heated_channel cases
#
# References:
#
# Moser, Kim, and Mansour. DNS of Turbulent Channel Flow up to Re_tau=590.
# Physics of Fluids, vol 11, 943-945, 1999.
#
# Kim, Moin & Moser. [Numerical Method, Re_tau = 178.12]. J. Fluid Mech.
# vol 177, 133-166, 1987.

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import fdsplotlib

plot_style = fdsplotlib.get_plot_style('fds')

outdir = '../../../out/Heated_Channel_Flow/'
expdir = '../../../exp/Heated_Channel_Flow/'
pltdir = '../../Manuals/FDS_Verification_Guide/SCRIPT_FIGURES/'

git_file = outdir + 'heated_channel_Pr_0p10_16_git.txt'
version_string = fdsplotlib.get_version_string(git_file)

# Plot DNS results

M = pd.read_csv(expdir + 'heated_channel_dns_data.csv', skiprows=0)

yp_up_mean = M['yp_up_mean'].values
up_mean = M['up_mean'].values
yp_Tp_mean = M['yp_Tp_mean'].values
Tp_mean_Pr0p10 = M['Tp_mean_Pr0p1'].values
Tp_mean_Pr0p71 = M['Tp_mean_Pr0p71'].values
Tp_mean_Pr2p00 = M['Tp_mean_Pr2p0'].values
range_indices = slice(0, len(yp_up_mean), 3)

fig1 = fdsplotlib.plot_to_fig(x_data=yp_up_mean, y_data=up_mean, marker_style='ro', data_label=r'DNS Re$_\tau=180$',
                              x_min=1, x_max=1000, y_min=0, y_max=30,
                              plot_type='semilogx',
                              revision_label=version_string,
                              x_label='$z^+$',
                              y_label='$u^+$')

fig2 = fdsplotlib.plot_to_fig(x_data=yp_Tp_mean[range_indices], y_data=Tp_mean_Pr0p10[range_indices], marker_style='bo', data_label=r'DNS Re$_\tau=180$',
                              x_min=1, x_max=1000, y_min=0, y_max=30,
                              plot_type='semilogx',
                              revision_label=version_string,
                              x_label='$z^+$',
                              y_label='$T^+$')

fdsplotlib.plot_to_fig(x_data=yp_Tp_mean[range_indices], y_data=Tp_mean_Pr0p71[range_indices], marker_style='bo', figure_handle=fig2)
fdsplotlib.plot_to_fig(x_data=yp_Tp_mean[range_indices], y_data=Tp_mean_Pr2p00[range_indices], marker_style='bo', figure_handle=fig2)

a2 = fig2.gca()
a2.text(225, 26, 'Pr=2.0',  fontsize=plot_style['Label_Font_Size'])
a2.text(225, 15, 'Pr=0.71', fontsize=plot_style['Label_Font_Size'])
a2.text(225, 5,  'Pr=0.10', fontsize=plot_style['Label_Font_Size'])

# plot the FDS results

H = 2
dpdx = 9.0088e-6
tau_w = 0.5 * dpdx * H
cp = 1
mu = 1.8216e-5
T_w = 20

devcfile = ['heated_channel_Pr_0p10_16_devc.csv',
            'heated_channel_Pr_0p71_16_devc.csv',
            'heated_channel_Pr_1p00_16_devc.csv',
            'heated_channel_Pr_2p00_16_devc.csv']

linefile = ['heated_channel_Pr_0p10_16_line.csv',
            'heated_channel_Pr_0p71_16_line.csv',
            'heated_channel_Pr_1p00_16_line.csv',
            'heated_channel_Pr_2p00_16_line.csv']

err = np.zeros(4) 

for i in range(4):

    M = pd.read_csv(outdir + devcfile[i], skiprows=1)

    rho = M['RHO'].iloc[-1]  # should be about 1.19
    u_tau = np.sqrt(tau_w / rho)
    delta_nu = (mu / rho) / u_tau

    # Find columns starting with HF0B and HF0T
    j1 = [j for j, col in enumerate(M.columns) if col == 'HF0B'][0]
    j2 = [j for j, col in enumerate(M.columns) if col == 'HF0T'][0]
    q_w = M.iloc[-1, j1:j2+1].mean()
    T_tau = q_w / (rho * u_tau * cp)

    M = pd.read_csv(outdir + linefile[i], skiprows=1)

    zp = M.iloc[0:16, 0].values / delta_nu

    if i in [0, 1, 3]:
        j1 = [j for j, col in enumerate(M.columns) if col == 'T11'][0]
        j2 = [j for j, col in enumerate(M.columns) if col == 'T115'][0]
        T1 = M.iloc[0:16, j1:j2+1].mean(axis=1).values      # bottom wall
        T2 = M.iloc[31:15:-1, j1:j2+1].mean(axis=1).values  # top wall
        Tp = (0.5 * (T1 + T2) - T_w) / T_tau
        if i==0:
            fdsplotlib.plot_to_fig(x_data=zp, y_data=Tp, marker_style='ksq-', data_label='FDS', figure_handle=fig2)
        else:
            fdsplotlib.plot_to_fig(x_data=zp, y_data=Tp, marker_style='ksq-', figure_handle=fig2)
    elif i == 2:
        j1 = [j for j, col in enumerate(M.columns) if col == 'U11'][0]
        j2 = [j for j, col in enumerate(M.columns) if col == 'U115'][0]
        u1 = M.iloc[0:16, j1:j2+1].mean(axis=1).values     # bottom wall
        u2 = M.iloc[31:15:-1, j1:j2+1].mean(axis=1).values # top wall
        up = 0.5 * (u1 + u2) / u_tau
        fdsplotlib.plot_to_fig(x_data=zp, y_data=up, marker_style='ksq-', data_label='FDS', figure_handle=fig1)

    if i == 0:
        err[i] = abs(np.mean(Tp) - np.nanmean(Tp_mean_Pr0p10)) / np.nanmean(Tp_mean_Pr0p10)
    elif i == 1:
        err[i] = abs(np.mean(Tp) - np.nanmean(Tp_mean_Pr0p71)) / np.nanmean(Tp_mean_Pr0p71)
    elif i == 2:
        err[i] = abs(np.mean(up) - np.nanmean(up_mean)) / np.nanmean(up_mean)
    elif i == 3:
        err[i] = abs(np.mean(Tp) - np.nanmean(Tp_mean_Pr2p00)) / np.nanmean(Tp_mean_Pr2p00)

    if err[i] > 1:
        print('Python Warning: heated_channel case ' + str(i+1) + ' out of tolerance')

fig1.savefig(pltdir + 'heated_channel_uplus.pdf', format='pdf')
fig2.savefig(pltdir + 'heated_channel_Tplus.pdf', format='pdf')

