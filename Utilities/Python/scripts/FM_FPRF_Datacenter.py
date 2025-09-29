
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import fdsplotlib

plot_style = fdsplotlib.get_plot_style('fds')

expdir = '../../../exp/FM_FPRF_Datacenter/'
outdir = '../../../out/FM_FPRF_Datacenter/'
pltdir = '../../Manuals/FDS_Validation_Guide/SCRIPT_FIGURES/FM_FPRF_Datacenter/'

git_file = outdir + 'FM_Datacenter_Veltest_High_git.txt'
version_string = fdsplotlib.get_version_string(git_file)

lb = [-2.0,-5.0,-1.0, 0.0,-0.6,-1.5,-0.4, 0.0]
ub = [ 2.0, 5.0, 3.0, 5.0, 0.6, 1.5, 0.8, 1.4]
dir1 = ['U-','V-','W-','Total ','U-','V-','W-','Total ']
dir2 = ['u','v','w','vel','u','v','w','vel']
exp_index = [0,1,2,6,0,1,2,6]
exp_rms_index = [3,4,5,7,3,4,5,7]

# Error lines for plots

x_err = np.arange(-10, 11, 1)
y_err = (0.1773**2 + (0.05*x_err)**2 + (0.06*x_err)**2)**0.5
y_err_p = x_err + 2*y_err
y_err_m = x_err - 2*y_err

# Make 8 plots, 4 for high flow and 4 for low flow

for j in range(8):

    if j == 0:  # Compute average velocity for high flow experiments

        exp_data = np.loadtxt(expdir + 'fm_datacenter_veltest_high.csv', delimiter=',', skiprows=1)
        fds_data = np.loadtxt(outdir + 'FM_Datacenter_Veltest_High_devc.csv', delimiter=',', skiprows=13)
        n_fds_data = fds_data.shape[0]
        fds_avg = np.zeros(195)
        for i in range(13, 208):
            fds_avg[i-13] = np.mean(fds_data[:, i-1])
        fds_u = fds_avg[0:65]
        fds_v = fds_avg[65:130]
        fds_w = fds_avg[130:195]
        fds_u_rms = fds_data[n_fds_data-1, 207:272]
        fds_v_rms = fds_data[n_fds_data-1, 272:337]
        fds_w_rms = fds_data[n_fds_data-1, 337:402]
        fds_tot = (fds_u**2 + fds_v**2 + fds_w**2)**0.5
        fds_tot_rms = (((fds_u * fds_u_rms)**2 + (fds_v * fds_v_rms)**2 + (fds_w * fds_w_rms)**2) / fds_tot**2)**0.5
        flow = 'High'
        title = 'High Flow Case'
        
    elif j==4:  # Compute average velocity for low flow experiments

        exp_data = np.loadtxt(expdir + 'fm_datacenter_veltest_low.csv', delimiter=',', skiprows=1)
        fds_data = np.loadtxt(outdir + 'FM_Datacenter_Veltest_Low_devc.csv', delimiter=',', skiprows=13)
        n_fds_data = fds_data.shape[0]
        fds_avg = np.zeros(183)
        for i in range(13, 196):
            fds_avg[i-13] = np.mean(fds_data[:, i-1])
        fds_u = fds_avg[0:61]
        fds_v = fds_avg[61:122]
        fds_w = fds_avg[122:183]
        fds_u_rms = fds_data[n_fds_data-1, 195:256]
        fds_v_rms = fds_data[n_fds_data-1, 256:317]
        fds_w_rms = fds_data[n_fds_data-1, 317:378]
        fds_tot = (fds_u**2 + fds_v**2 + fds_w**2)**0.5
        fds_tot_rms = (((fds_u * fds_u_rms)**2 + (fds_v * fds_v_rms)**2 + (fds_w * fds_w_rms)**2) / fds_tot**2)**0.5
        flow='Low'
        title = 'Low Flow Case'
        
    # Choose the proper FDS values
    
    if j==0 or j==4:
        fds_points = fds_u
        fds_error  = fds_u_rms
    elif j==1 or j==5:
        fds_points = fds_v
        fds_error  = fds_v_rms
    elif j==2 or j==6:
        fds_points = fds_w
        fds_error  = fds_w_rms
    elif j==3 or j==7:
        fds_points = fds_tot
        fds_error  = fds_tot_rms

    # Plot and save the figure

    fig = fdsplotlib.plot_to_fig(x_data=x_err, y_data=x_err,
                                 marker_style='k-',
                                 x_min=lb[j], x_max=ub[j], y_min=lb[j], y_max=ub[j],
                                 figure_size=(plot_style['Scat_Paper_Width'],plot_style['Scat_Paper_Height']),
                                 plot_size=(plot_style['Scat_Plot_Width'],plot_style['Scat_Plot_Height']),
                                 plot_origin=(plot_style['Scat_Plot_X'],plot_style['Scat_Plot_Y']),
                                 revision_label=version_string, plot_title=title,
                                 x_label=f'Measured {dir1[j]}Velocity (m/s)',
                                 y_label=f'Predicted {dir1[j]}Velocity (m/s)')
    fdsplotlib.plot_to_fig(x_data=x_err, y_data=y_err_p, marker_style='k--', figure_handle=fig)
    fdsplotlib.plot_to_fig(x_data=x_err, y_data=y_err_m, marker_style='k--', figure_handle=fig)
    fdsplotlib.plot_to_fig(x_data=exp_data[:,exp_index[j]], y_data=fds_points,
                           x_error=exp_data[:, exp_rms_index[j]], y_error=fds_u_rms,
                           marker_style='ro', figure_handle=fig)

    plt.savefig(pltdir + f'FM_Datacenter_Veltest_{flow}_{dir2[j]}.pdf', format='pdf')
    plt.close()

