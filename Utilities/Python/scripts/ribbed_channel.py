
# Casara Art Ribbed Channel cases

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import sys
import fdsplotlib

plot_style = fdsplotlib.get_plot_style('fds')

outdir = '../../../out/Casara_Arts_Ribbed_Channel/'
expdir = '../../../exp/Casara_Arts_Ribbed_Channel/'
pltdir = '../../Manuals/FDS_Verification_Guide/SCRIPT_FIGURES/'

git_file = outdir + 'ribbed_channel_20_git.txt'
version_string = fdsplotlib.get_version_string(git_file)

nx = [20, 40, 80, 160]
lnx = len(nx)
Ub = 6.2  # exp bulk velocity
L = 0.3   # channel length
D = 0.1   # channel height
h = 0.03
dx = L / np.array([10, 20, 40, 80])
fds_marker = ['r+-', 'c^-', 'g>-', 'k-']
fds_key = [r'FDS $h/\delta x=3$', r'FDS $h/\delta x=6$', r'FDS $h/\delta x=12$', r'FDS $h/\delta x=24$']
geom = ['_', '_geom_']

DATA = pd.read_csv(os.path.join(expdir, 'ribbed_channel_data.csv'))

# Main loop over regular and complex geometry types

for ii in range(len(geom)):
    
    # Bulk velocity 

    M = {}
    for i in range(lnx):
        filename = os.path.join(outdir, f'ribbed_channel{geom[ii]}{nx[i]}_devc.csv')
        M[i] = pd.read_csv(filename, skiprows=1)
    
    fig = fdsplotlib.plot_to_fig(x_data=[0,100], y_data=[Ub,Ub], marker_style='k-', data_label='Exp Bulk Velocity',
                                 x_min=0, x_max=10, y_min=0, y_max=1.5*Ub,
                                 revision_label=version_string,
                                 x_label='Time (s)',
                                 y_label='Bulk Velocity (m/s)')

    for i in range(lnx):
        t_fds = M[i].iloc[:, 0].values
        Ub_fds = M[i].iloc[:, 1].values
        fdsplotlib.plot_to_fig(x_data=t_fds, y_data=Ub_fds, figure_handle=fig, marker_style=fds_marker[i], data_label=fds_key[i])
        t_range = np.where(t_fds > 2)[0]
        if abs(np.mean(Ub_fds[t_range]) - Ub) / Ub > 0.01:
            print(f'Matlab Warning: Ub mean nx {geom[ii]}{nx[i]} = {np.mean(Ub_fds[t_range])}')
    
    plt.savefig(os.path.join(pltdir, f'ribbed_channel{geom[ii]}Ubulk.pdf'), format='pdf')
    plt.close()
    
    # Streamwise U along bottom of channel

    for i in range(lnx):
        filename = os.path.join(outdir, f'ribbed_channel{geom[ii]}{nx[i]}_line.csv')
        M[i] = pd.read_csv(filename, skiprows=1)
    
    j = DATA.columns.get_loc('x/h U strm')
    xoh = DATA.iloc[:, j].values
    j = DATA.columns.get_loc('U strm')
    u_data = DATA.iloc[:, j].values
    
    fig = fdsplotlib.plot_to_fig(x_data=xoh, y_data=u_data, marker_style='bo', data_label='PIV Data',
                                 x_min=0, x_max=10, y_min=-0.6, y_max=0.6,
                                 revision_label=version_string,
                                 x_label='$x/h$',
                                 y_label=r'$U/U_{\rm b}$')
    
    for i in range(lnx):
        j = M[i].columns.get_loc('u_strm_bot-x')
        x = M[i].iloc[:, j].values
        I = np.where(x < 0)[0]
        x[I] = x[I] + L
        I = np.argsort(x)
        x = x[I]
        j = M[i].columns.get_loc('u_strm_bot')
        u_fds = M[i].iloc[:, j].values
        u_fds = u_fds[I]
        fdsplotlib.plot_to_fig(x_data=x/h, y_data=u_fds/Ub, figure_handle=fig, marker_style=fds_marker[i], data_label=fds_key[i])
    
    plt.savefig(os.path.join(pltdir, f'ribbed_channel{geom[ii]}u_strm.pdf'), format='pdf')
    plt.close()
    
    # Streamwise urms along bottom of channel

    j = DATA.columns.get_loc('x/h urms strm')
    xoh = DATA.iloc[:, j].values
    j = DATA.columns.get_loc('urms strm')
    urms_data = DATA.iloc[:, j].values
    
    fig = fdsplotlib.plot_to_fig(x_data=xoh, y_data=urms_data, marker_style='bo', data_label='PIV Data',
                                 x_min=0, x_max=10, y_min=0, y_max=0.6,
                                 revision_label=version_string,
                                 x_label='$x/h$',
                                 y_label=r'$U_{\rm rms}/U_{\rm b}$')
    
    for i in range(lnx):
        j = M[i].columns.get_loc('urms_strm_bot-x')
        x = M[i].iloc[:, j].values
        I = np.where(x < 0)[0]
        x[I] = x[I] + L
        I = np.argsort(x)
        x = x[I]
        j = M[i].columns.get_loc('urms_strm_bot')
        urms_fds = M[i].iloc[:, j].values
        urms_fds = urms_fds[I]
        fdsplotlib.plot_to_fig(x_data=x/h, y_data=urms_fds/Ub, figure_handle=fig, marker_style=fds_marker[i], data_label=fds_key[i])
    
    plt.savefig(os.path.join(pltdir, f'ribbed_channel{geom[ii]}urms_strm.pdf'), format='pdf')
    plt.close()
    
    # Streamwise U profile at x/h=0 (center of rib)

    j = DATA.columns.get_loc('y/h U prof')
    yoh = DATA.iloc[:, j].values
    j = DATA.columns.get_loc('U prof')
    u_data = DATA.iloc[:, j].values
    
    fig = fdsplotlib.plot_to_fig(x_data=u_data, y_data=yoh, marker_style='bo', data_label='PIV Data',
                                 x_min=0, x_max=2, y_min=1, y_max=D/h,
                                 revision_label=version_string,
                                 x_label='$z/h$',
                                 y_label=r'$U/U_{\rm b}$')
    
    for i in range(lnx):
        j = M[i].columns.get_loc('u_prof_rib-z')
        I = np.where((M[i].iloc[:, j].values > h) & (M[i].iloc[:, j].values < D))[0]
        y = M[i].iloc[I, j].values
        j = M[i].columns.get_loc('u_prof_rib')
        u_fds = M[i].iloc[I, j].values
        fdsplotlib.plot_to_fig(x_data=u_fds/Ub, y_data=y/h, figure_handle=fig, marker_style=fds_marker[i], data_label=fds_key[i])
    
    plt.savefig(os.path.join(pltdir, f'ribbed_channel{geom[ii]}u_prof.pdf'), format='pdf')
    plt.close()
    
    # Streamwise urms profile at x/h=0 (center of rib)

    j = DATA.columns.get_loc('y/h urms prof')
    yoh = DATA.iloc[:, j].values
    j = DATA.columns.get_loc('urms prof')
    urms_data = DATA.iloc[:, j].values
    
    fig = fdsplotlib.plot_to_fig(x_data=urms_data, y_data=yoh, marker_style='bo', data_label='PIV Data',
                                 x_min=0, x_max=1, y_min=1, y_max=D/h,
                                 revision_label=version_string,
                                 x_label='$z/h$',
                                 y_label=r'$U_{\rm rms}/U_{\rm b}$')
    
    for i in range(lnx):
        j = M[i].columns.get_loc('urms_prof_rib-z')
        I = np.where((M[i].iloc[:, j].values > h) & (M[i].iloc[:, j].values < D))[0]
        y = M[i].iloc[I, j].values
        j = M[i].columns.get_loc('urms_prof_rib')
        urms_fds = M[i].iloc[I, j].values
        fdsplotlib.plot_to_fig(x_data=urms_fds/Ub, y_data=y/h, figure_handle=fig, marker_style=fds_marker[i], data_label=fds_key[i])
    
    plt.savefig(os.path.join(pltdir, f'ribbed_channel{geom[ii]}urms_prof.pdf'), format='pdf')
    plt.close()

