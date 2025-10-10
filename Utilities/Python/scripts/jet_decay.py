
# Verification Guide, Jet Centerline Velocity Decay

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import rc
import os
import fdsplotlib

plot_style = fdsplotlib.get_plot_style('fds')

outdir = '../../../out/Turbulent_Jet/'
pltdir = '../../Manuals/FDS_Verification_Guide/SCRIPT_FIGURES/'

git_file = outdir + 'jet_csmag_dx10cm_git.txt'
version_string = fdsplotlib.get_version_string(git_file)

# gather FDS results

M = pd.read_csv(outdir + 'jet_csmag_dx10cm_line.csv', skiprows=2)
u_csmag_10 = M.iloc[:, 1].values

M = pd.read_csv(outdir + 'jet_dsmag_dx10cm_line.csv', skiprows=2)
u_dsmag_10 = M.iloc[:, 1].values

M = pd.read_csv(outdir + 'jet_deardorff_dx10cm_line.csv', skiprows=2)
u_deardorff_10 = M.iloc[:, 1].values

M = pd.read_csv(outdir + 'jet_vreman_dx10cm_line.csv', skiprows=2)
u_vreman_10 = M.iloc[:, 1].values

M = pd.read_csv(outdir + 'jet_csmag_dx5cm_line.csv', skiprows=2)
u_csmag_5 = M.iloc[:, 1].values

M = pd.read_csv(outdir + 'jet_dsmag_dx5cm_line.csv', skiprows=2)
u_dsmag_5 = M.iloc[:, 1].values

M = pd.read_csv(outdir + 'jet_deardorff_dx5cm_line.csv', skiprows=2)
u_deardorff_5 = M.iloc[:, 1].values

M = pd.read_csv(outdir + 'jet_vreman_dx5cm_line.csv', skiprows=2)
u_vreman_5 = M.iloc[:, 1].values

# analytical solutions

x = M.iloc[:, 0].values
u_0 = 2.1
h = 0.8
b = 0.8
m_1 = 0.12
m_2 = 0.20

u_1 = np.zeros(len(x))
for j in range(len(x)):
    if x[j] < h/m_1:
        u_1[j] = u_0
    else:
        u_1[j] = u_0/(m_1*x[j])*np.sqrt(b*h)

u_2 = np.zeros(len(x))
for j in range(len(x)):
    if x[j] < h/m_2:
        u_2[j] = u_0
    else:
        u_2[j] = u_0/(m_2*x[j])*np.sqrt(b*h)

# Create figure

fig = fdsplotlib.plot_to_fig(x_data=[-1,-1], y_data=[-1,-1],
                             x_min=0, x_max=25, y_min=0, y_max=1.2,
                             revision_label=version_string,
                             plot_title='Jet Centerline Velocity Decay',
                             legend_location='upper right',
                             legend_fontsize=8,
                             x_label='$x/h$',
                             y_label='$u_{\\hbox{\\tiny max}}/u_0$')

fdsplotlib.plot_to_fig(x_data=x/h, y_data=u_1/u_0, figure_handle=fig, marker_style='k--', data_label='analytical, $m=0.12$')
fdsplotlib.plot_to_fig(x_data=x/h, y_data=u_2/u_0, figure_handle=fig, marker_style='k-' , data_label='analytical, $m=0.20$')
fdsplotlib.plot_to_fig(x_data=x/h, y_data=u_csmag_10/u_0, figure_handle=fig, marker_style='g--', data_label='csmag, $h/\\delta x=8$')
fdsplotlib.plot_to_fig(x_data=x/h, y_data=u_csmag_5 /u_0, figure_handle=fig, marker_style='g-' , data_label='csmag, $h/\\delta x=16$')
fdsplotlib.plot_to_fig(x_data=x/h, y_data=u_dsmag_10/u_0, figure_handle=fig, marker_style='c--', data_label='dsmag, $h/\\delta x=8$')
fdsplotlib.plot_to_fig(x_data=x/h, y_data=u_dsmag_5 /u_0, figure_handle=fig, marker_style='c-' , data_label='dsmag, $h/\\delta x=16$')
fdsplotlib.plot_to_fig(x_data=x/h, y_data=u_deardorff_10/u_0, figure_handle=fig, marker_style='b--', data_label='Deardorff, $h/\\delta x=8$')
fdsplotlib.plot_to_fig(x_data=x/h, y_data=u_deardorff_5 /u_0, figure_handle=fig, marker_style='b-' , data_label='Deardorff, $h/\\delta x=16$')
fdsplotlib.plot_to_fig(x_data=x/h, y_data=u_vreman_10/u_0, figure_handle=fig, marker_style='r--', data_label='Vreman, $h/\\delta x=8$')
fdsplotlib.plot_to_fig(x_data=x/h, y_data=u_vreman_5 /u_0, figure_handle=fig, marker_style='r-' , data_label='Vreman, $h/\\delta x=16$')

plt.savefig(pltdir + 'jet_decay.pdf', format='pdf')
plt.close()

