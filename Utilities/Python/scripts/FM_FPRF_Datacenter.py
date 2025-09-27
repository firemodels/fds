
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import warnings
import os
from matplotlib.backends.backend_pdf import PdfPages

# include FDS plot styles, etc.
import fdsplotlib

# Get plot style parameters
plot_style = fdsplotlib.get_plot_style('fds')

expdir = '../../../exp/FM_FPRF_Datacenter/'
outdir = '../../../out/FM_FPRF_Datacenter/'

git_file = outdir + 'FM_Datacenter_Veltest_High_git.txt'
version_string = fdsplotlib.get_version_string(git_file)

# High flow test
exp_data = np.loadtxt(expdir + 'fm_datacenter_veltest_high.csv', delimiter=',', skiprows=1)
fds_data = np.loadtxt(outdir + 'FM_Datacenter_Veltest_High_devc.csv', delimiter=',', skiprows=13)
n_fds_data = fds_data.shape[0]

# compute average velocity
fds_avg = np.zeros(195)  # 207-12 = 195
for i in range(13, 208):  # MATLAB 13:207 -> Python 12:207 (inclusive)
    fds_avg[i-13] = np.mean(fds_data[:, i-1])  # MATLAB is 1-indexed, Python is 0-indexed

fds_u = fds_avg[0:65]  # MATLAB 1:65 -> Python 0:65
fds_v = fds_avg[65:130]  # MATLAB 66:130 -> Python 65:130
fds_w = fds_avg[130:195]  # MATLAB 131:195 -> Python 130:195
fds_u_rms = fds_data[n_fds_data-1, 195:260]  # MATLAB 196:260 -> Python 195:260
fds_v_rms = fds_data[n_fds_data-1, 260:325]  # MATLAB 261:325 -> Python 260:325
fds_w_rms = fds_data[n_fds_data-1, 325:390]  # MATLAB 326:390 -> Python 325:390
fds_tot = (fds_u**2 + fds_v**2 + fds_w**2)**0.5
fds_tot_rms = (((fds_u * fds_u_rms)**2 + (fds_v * fds_v_rms)**2 + (fds_w * fds_w_rms)**2) / fds_tot**2)**0.5

# set plot error lines
x_err = np.arange(-10, 11, 1)  # MATLAB -10:1:10
y_err = (0.1773**2 + (0.05*x_err)**2 + (0.06*x_err)**2)**0.5
y_err_p = x_err + 2*y_err
y_err_m = x_err - 2*y_err

# Create figure and plot U-velocity comparison

fig = fdsplotlib.plot_to_fig(x_data=x_err, y_data=x_err,
                             marker_style='k-',
                             x_min=-2, x_max=2, y_min=-2, y_max=2,
                             figure_size=(plot_style['Scat_Paper_Width'],plot_style['Scat_Paper_Height']),
                             plot_size=(plot_style['Scat_Plot_Width'],plot_style['Scat_Plot_Height']),
                             plot_origin=(plot_style['Scat_Plot_X'],plot_style['Scat_Plot_Y']),
                             revision_label=version_string,
                             plot_title='High Flow',
                             x_label='Measured U-Velocity (m/s)',
                             y_label='Predicted U-Velocity (m/s)')
fdsplotlib.plot_to_fig(x_data=x_err, y_data=y_err_p, marker_style='k--', figure_handle=fig)
fdsplotlib.plot_to_fig(x_data=x_err, y_data=y_err_m, marker_style='k--', figure_handle=fig)
fdsplotlib.plot_to_fig(x_data=exp_data[:, 0], y_data=fds_u, marker_style='ro', figure_handle=fig)
#herrorbar(exp_data[:, 0], fds_u, exp_data[:, 3], 'ro')
#plt.errorbar(exp_data[:, 0], fds_u, yerr=fds_u_rms, fmt='ro')
plotname = '../../Manuals/FDS_Validation_Guide/SCRIPT_FIGURES/FM_FPRF_Datacenter/FM_Datacenter_Veltest_High_u'
plt.savefig(plotname + '.pdf', format='pdf')
plt.close()

# Create figure and plot V-velocity comparison

fig = fdsplotlib.plot_to_fig(x_data=x_err, y_data=x_err,
                             marker_style='k-',
                             x_min=-5, x_max=5, y_min=-5, y_max=5,
                             figure_size=(plot_style['Scat_Paper_Width'],plot_style['Scat_Paper_Height']),
                             plot_size=(plot_style['Scat_Plot_Width'],plot_style['Scat_Plot_Height']),
                             plot_origin=(plot_style['Scat_Plot_X'],plot_style['Scat_Plot_Y']),
                             revision_label=version_string,
                             plot_title='High Flow',
                             x_label='Measured V-Velocity (m/s)',
                             y_label='Predicted V-Velocity (m/s)')
fdsplotlib.plot_to_fig(x_data=x_err, y_data=y_err_p, marker_style='k--', figure_handle=fig)
fdsplotlib.plot_to_fig(x_data=x_err, y_data=y_err_m, marker_style='k--', figure_handle=fig)
fdsplotlib.plot_to_fig(x_data=exp_data[:, 1], y_data=fds_v, marker_style='ro', figure_handle=fig)
#herrorbar(exp_data[:, 1], fds_v, exp_data[:, 4], 'ro')
#plt.errorbar(exp_data[:, 1], fds_v, yerr=fds_v_rms, fmt='ro')
plotname = '../../Manuals/FDS_Validation_Guide/SCRIPT_FIGURES/FM_FPRF_Datacenter/FM_Datacenter_Veltest_High_v'
plt.savefig(plotname + '.pdf', format='pdf')
plt.close()

# Create figure and plot W-velocity comparison

fig = fdsplotlib.plot_to_fig(x_data=x_err, y_data=x_err,
                             marker_style='k-',
                             x_min=-1, x_max=3, y_min=-1, y_max=3,
                             figure_size=(plot_style['Scat_Paper_Width'],plot_style['Scat_Paper_Height']),
                             plot_size=(plot_style['Scat_Plot_Width'],plot_style['Scat_Plot_Height']),
                             plot_origin=(plot_style['Scat_Plot_X'],plot_style['Scat_Plot_Y']),
                             revision_label=version_string,
                             plot_title='High Flow',
                             x_label='Measured W-Velocity (m/s)',
                             y_label='Predicted W-Velocity (m/s)')
fdsplotlib.plot_to_fig(x_data=x_err, y_data=y_err_p, marker_style='k--', figure_handle=fig)
fdsplotlib.plot_to_fig(x_data=x_err, y_data=y_err_m, marker_style='k--', figure_handle=fig)
fdsplotlib.plot_to_fig(x_data=exp_data[:, 2], y_data=fds_w, marker_style='ro', figure_handle=fig)
#herrorbar(exp_data[:, 2], fds_w, exp_data[:, 5], 'ro')
#plt.errorbar(exp_data[:, 2], fds_w, yerr=fds_w_rms, fmt='ro')
plotname = '../../Manuals/FDS_Validation_Guide/SCRIPT_FIGURES/FM_FPRF_Datacenter/FM_Datacenter_Veltest_High_w'
plt.savefig(plotname + '.pdf', format='pdf')
plt.close()

# Create figure and plot total velocity comparison

fig = fdsplotlib.plot_to_fig(x_data=x_err, y_data=x_err,
                             marker_style='k-',
                             x_min=0, x_max=5, y_min=0, y_max=5,
                             figure_size=(plot_style['Scat_Paper_Width'],plot_style['Scat_Paper_Height']),
                             plot_size=(plot_style['Scat_Plot_Width'],plot_style['Scat_Plot_Height']),
                             plot_origin=(plot_style['Scat_Plot_X'],plot_style['Scat_Plot_Y']),
                             revision_label=version_string,
                             plot_title='High Flow',
                             x_label='Measured Total Velocity (m/s)',
                             y_label='Predicted Total Velocity (m/s)')
fdsplotlib.plot_to_fig(x_data=x_err, y_data=y_err_p, marker_style='k--', figure_handle=fig)
fdsplotlib.plot_to_fig(x_data=x_err, y_data=y_err_m, marker_style='k--', figure_handle=fig)
fdsplotlib.plot_to_fig(x_data=exp_data[:, 6], y_data=fds_tot, marker_style='ro', figure_handle=fig)
#herrorbar(exp_data[:, 6], fds_tot, exp_data[:, 7], 'ro')
#plt.errorbar(exp_data[:, 6], fds_tot, yerr=fds_tot_rms, fmt='ro')
plotname = '../../Manuals/FDS_Validation_Guide/SCRIPT_FIGURES/FM_FPRF_Datacenter/FM_Datacenter_Veltest_High_vel'
plt.savefig(plotname + '.pdf', format='pdf')
plt.close()


# Low Flow Tests

exp_data = np.loadtxt(expdir + 'fm_datacenter_veltest_low.csv', delimiter=',', skiprows=1)
fds_data = np.loadtxt(outdir + 'FM_Datacenter_Veltest_Low_devc.csv', delimiter=',', skiprows=13)
n_fds_data = fds_data.shape[0]

# compute average velocity
fds_avg = np.zeros(183)  # 195-12 = 183
for i in range(13, 196):  # MATLAB 13:195 -> Python 12:195 (inclusive)
    fds_avg[i-13] = np.mean(fds_data[:, i-1])  # MATLAB is 1-indexed, Python is 0-indexed

fds_u = fds_avg[0:61]  # MATLAB 1:61 -> Python 0:61
fds_v = fds_avg[61:122]  # MATLAB 62:122 -> Python 61:122
fds_w = fds_avg[122:183]  # MATLAB 123:183 -> Python 122:183
fds_u_rms = fds_data[n_fds_data-1, 183:244]  # MATLAB 184:244 -> Python 183:244
fds_v_rms = fds_data[n_fds_data-1, 244:305]  # MATLAB 245:305 -> Python 244:305
fds_w_rms = fds_data[n_fds_data-1, 305:366]  # MATLAB 306:366 -> Python 305:366
fds_tot = (fds_u**2 + fds_v**2 + fds_w**2)**0.5
fds_tot_rms = (((fds_u * fds_u_rms)**2 + (fds_v * fds_v_rms)**2 + (fds_w * fds_w_rms)**2) / fds_tot**2)**0.5

# Create figure and plot U-velocity comparison, low flow

fig = fdsplotlib.plot_to_fig(x_data=x_err, y_data=x_err,
                             marker_style='k-',
                             x_min=-0.6, x_max=0.6, y_min=-0.6, y_max=0.6,
                             figure_size=(plot_style['Scat_Paper_Width'],plot_style['Scat_Paper_Height']),
                             plot_size=(plot_style['Scat_Plot_Width'],plot_style['Scat_Plot_Height']),
                             plot_origin=(plot_style['Scat_Plot_X'],plot_style['Scat_Plot_Y']),
                             revision_label=version_string,
                             plot_title='Low Flow',
                             x_label='Measured U-Velocity (m/s)',
                             y_label='Predicted U-Velocity (m/s)')
fdsplotlib.plot_to_fig(x_data=x_err, y_data=y_err_p, marker_style='k--', figure_handle=fig)
fdsplotlib.plot_to_fig(x_data=x_err, y_data=y_err_m, marker_style='k--', figure_handle=fig)
fdsplotlib.plot_to_fig(x_data=exp_data[:, 0], y_data=fds_u, marker_style='ro', figure_handle=fig)
#herrorbar(exp_data[:, 0], fds_u, exp_data[:, 3], 'ro')
#plt.errorbar(exp_data[:, 0], fds_u, yerr=fds_u_rms, fmt='ro')
plotname = '../../Manuals/FDS_Validation_Guide/SCRIPT_FIGURES/FM_FPRF_Datacenter/FM_Datacenter_Veltest_Low_u'
plt.savefig(plotname + '.pdf', format='pdf')
plt.close()

# Create figure and plot V-velocity comparison, low flow

fig = fdsplotlib.plot_to_fig(x_data=x_err, y_data=x_err,
                             marker_style='k-',
                             x_min=-1.5, x_max=1.5, y_min=-1.5, y_max=1.5,
                             figure_size=(plot_style['Scat_Paper_Width'],plot_style['Scat_Paper_Height']),
                             plot_size=(plot_style['Scat_Plot_Width'],plot_style['Scat_Plot_Height']),
                             plot_origin=(plot_style['Scat_Plot_X'],plot_style['Scat_Plot_Y']),
                             revision_label=version_string,
                             plot_title='Low Flow',
                             x_label='Measured V-Velocity (m/s)',
                             y_label='Predicted V-Velocity (m/s)')
fdsplotlib.plot_to_fig(x_data=x_err, y_data=y_err_p, marker_style='k--', figure_handle=fig)
fdsplotlib.plot_to_fig(x_data=x_err, y_data=y_err_m, marker_style='k--', figure_handle=fig)
fdsplotlib.plot_to_fig(x_data=exp_data[:, 1], y_data=fds_v, marker_style='ro', figure_handle=fig)
#herrorbar(exp_data[:, 1], fds_v, exp_data[:, 4], 'ro')
#plt.errorbar(exp_data[:, 1], fds_v, yerr=fds_v_rms, fmt='ro')
plotname = '../../Manuals/FDS_Validation_Guide/SCRIPT_FIGURES/FM_FPRF_Datacenter/FM_Datacenter_Veltest_Low_v'
plt.savefig(plotname + '.pdf', format='pdf')
plt.close()

# Create figure and plot W-velocity comparison, low flow

fig = fdsplotlib.plot_to_fig(x_data=x_err, y_data=x_err,
                             marker_style='k-',
                             x_min=-0.4, x_max=0.8, y_min=-0.4, y_max=0.8,
                             figure_size=(plot_style['Scat_Paper_Width'],plot_style['Scat_Paper_Height']),
                             plot_size=(plot_style['Scat_Plot_Width'],plot_style['Scat_Plot_Height']),
                             plot_origin=(plot_style['Scat_Plot_X'],plot_style['Scat_Plot_Y']),
                             revision_label=version_string,
                             plot_title='Low Flow',
                             x_label='Measured W-Velocity (m/s)',
                             y_label='Predicted W-Velocity (m/s)')
fdsplotlib.plot_to_fig(x_data=x_err, y_data=y_err_p, marker_style='k--', figure_handle=fig)
fdsplotlib.plot_to_fig(x_data=x_err, y_data=y_err_m, marker_style='k--', figure_handle=fig)
fdsplotlib.plot_to_fig(x_data=exp_data[:, 2], y_data=fds_w, marker_style='ro', figure_handle=fig)
#herrorbar(exp_data[:, 2], fds_w, exp_data[:, 5], 'ro')
#plt.errorbar(exp_data[:, 2], fds_w, yerr=fds_w_rms, fmt='ro')
plotname = '../../Manuals/FDS_Validation_Guide/SCRIPT_FIGURES/FM_FPRF_Datacenter/FM_Datacenter_Veltest_Low_w'
plt.savefig(plotname + '.pdf', format='pdf')
plt.close()

# Create figure and plot total velocity comparison

fig = fdsplotlib.plot_to_fig(x_data=x_err, y_data=x_err,
                             marker_style='k-',
                             x_min=0, x_max=1.4, y_min=0, y_max=1.4,
                             figure_size=(plot_style['Scat_Paper_Width'],plot_style['Scat_Paper_Height']),
                             plot_size=(plot_style['Scat_Plot_Width'],plot_style['Scat_Plot_Height']),
                             plot_origin=(plot_style['Scat_Plot_X'],plot_style['Scat_Plot_Y']),
                             revision_label=version_string,
                             plot_title='Low Flow',
                             x_label='Measured Total Velocity (m/s)',
                             y_label='Predicted Total Velocity (m/s)')
fdsplotlib.plot_to_fig(x_data=x_err, y_data=y_err_p, marker_style='k--', figure_handle=fig)
fdsplotlib.plot_to_fig(x_data=x_err, y_data=y_err_m, marker_style='k--', figure_handle=fig)
fdsplotlib.plot_to_fig(x_data=exp_data[:, 6], y_data=fds_tot, marker_style='ro', figure_handle=fig)
#herrorbar(exp_data[:, 6], fds_tot, exp_data[:, 7], 'ro')
#plt.errorbar(exp_data[:, 6], fds_tot, yerr=fds_tot_rms, fmt='ro')
plotname = '../../Manuals/FDS_Validation_Guide/SCRIPT_FIGURES/FM_FPRF_Datacenter/FM_Datacenter_Veltest_Low_vel'
plt.savefig(plotname + '.pdf', format='pdf')
plt.close()


