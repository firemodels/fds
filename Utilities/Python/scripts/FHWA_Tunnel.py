
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
import os
import fdsplotlib

plot_style = fdsplotlib.get_plot_style('fds')

outdir = '../../../out/FHWA_Tunnel/'
expdir = '../../../exp/FHWA_Tunnel/'
pltdir = '../../Manuals/FDS_Validation_Guide/SCRIPT_FIGURES/FHWA_Tunnel/'

git_file = outdir + 'IFAB-07_cat_git.txt'
version_string = fdsplotlib.get_version_string(git_file)

pos = [1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.5, 9.5, 10.5, 11.5]
test = ['IFAB-07', 'IFAB-08', 'IFAB-09', 'IFAB-10', 'IFAB-11', 'IFAB-13', 'IFAB-14', 'IFAB-15', 'IFAB-19', 'IFAB-22', 'IFAB-24']
test2 = ['Test 7', 'Test 8', 'Test 9', 'Test 10', 'Test 11', 'Test 13', 'Test 14', 'Test 15', 'Test 19', 'Test 22', 'Test 24']
single_level = [50]
setpoint = [10000, 400, 399, 338, 240, 322, 390, 420, 360, 10000, 10000]

# For each experiment, make a contour plot of the extent of a single temperature contour at each time during the experiment

for k in range(11):  # Experiments

    fig = fdsplotlib.plot_to_fig(x_data=[5.5,5.5], y_data=[0,15], marker_style='k--',
                                 x_min=0, x_max=15, y_min=0, y_max=15,
                                 revision_label=version_string,
                                 plot_title=test2[k],
                                 x_label='Position (m)',
                                 y_label='Time (min)')
    fdsplotlib.plot_to_fig(x_data=[0.0,15.0], y_data=[setpoint[k]/60,setpoint[k]/60], figure_handle=fig, marker_style='k:')

    M_data = pd.read_csv(outdir + test[k] + '_cat_devc.csv', skiprows=2)
    M = {'data': M_data.values}
    E_data = pd.read_csv(expdir + test[k] + '_avg.csv', skiprows=2)
    E = {'data': E_data.values}

    X_mod, Y_mod = np.meshgrid(pos[0:16], M['data'][:, 0]/60)
    X_exp, Y_exp = np.meshgrid(pos[0:16], E['data'][:, 0]/60)
    
    Z_mod = np.zeros((len(M['data'][:, 0]), 16))
    Z_exp = np.zeros((len(E['data'][:, 0]), 16))
    
    for kk in range(len(M['data'][:, 0])):
        for ii in range(16):
            Z_mod[kk, ii] = M['data'][kk, ii+1]

    for kk in range(len(E['data'][:, 0])):
        for ii in range(16):
            Z_exp[kk, ii] = E['data'][kk, ii+4]

    newpoints = 100
    X_mod_interp, Y_mod_interp = np.meshgrid(
        np.linspace(np.min(X_mod), np.max(X_mod), newpoints),
        np.linspace(np.min(Y_mod), np.max(Y_mod), newpoints)
    )
    
    points_mod = np.column_stack((X_mod.ravel(), Y_mod.ravel()))
    values_mod = Z_mod.ravel()
    Z_mod_interp = griddata(points_mod, values_mod, (X_mod_interp, Y_mod_interp), method='cubic')

    X_exp_interp, Y_exp_interp = np.meshgrid(
        np.linspace(np.min(X_exp), np.max(X_exp), newpoints),
        np.linspace(np.min(Y_exp), np.max(Y_exp), newpoints)
    )
    
    points_exp = np.column_stack((X_exp.ravel(), Y_exp.ravel()))
    values_exp = Z_exp.ravel()
    Z_exp_interp = griddata(points_exp, values_exp, (X_exp_interp, Y_exp_interp), method='cubic')

    CS_mod = plt.contour(X_mod_interp, Y_mod_interp, Z_mod_interp, single_level, colors='red', linestyles='-')
    plt.clabel(CS_mod, fontsize=8, colors='red')
    
    CS_exp = plt.contour(X_exp_interp, Y_exp_interp, Z_exp_interp, single_level, colors='black', linestyles='-')
    plt.clabel(CS_exp, fontsize=8, colors='black')

    from matplotlib.lines import Line2D
    legend_elements = [ Line2D([0], [0], color='red', label='FDS'),
                        Line2D([0], [0], color='black', label='Exp') ]
    plt.legend(handles=legend_elements, loc='upper right', fontsize=plot_style['Key_Font_Size'])

    plt.savefig(pltdir + test[k] + '_tvT.pdf', format='pdf')

    plt.clf()  # Clear the figure for next iteration

