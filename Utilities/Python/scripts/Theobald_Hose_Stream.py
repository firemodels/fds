
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
import os
from pathlib import Path
import fdsplotlib

plot_style = fdsplotlib.get_plot_style('fds')

params_loc = '../../Validation/Theobald_Hose_Stream/FDS_Input_Files/Build_Input_Files/'
outdir = '../../../out/Theobald_Hose_Stream/'
expdir = '../../../exp/Theobald_Hose_Stream/'
plot_dir = '../../Manuals/FDS_Validation_Guide/SCRIPT_FIGURES/Theobald_Hose_Stream/'

git_file = outdir + 'Theobald_Test_0_git.txt'
version_string = fdsplotlib.get_version_string(git_file)

chid_loc = params_loc + 'paramfile.csv'
exp_loc = expdir + 'theobald_effect_1981_fds.csv'

# get list of FDS CHIDs
import_chid_data = pd.read_csv(chid_loc, header=0)
chid_list = import_chid_data.iloc[:, 1].tolist()
s_Theobald = len(chid_list)

# import exp results
Theobald_import_exp = pd.read_csv(exp_loc).values
max_range_exp = Theobald_import_exp[:, 7]  # Column 8 in MATLAB (0-indexed in Python)
max_height_exp = Theobald_import_exp[:, 5]  # Column 6 in MATLAB
max_height_dist_exp = Theobald_import_exp[:, 6]  # Column 7 in MATLAB

# import exp parameters
exp_noz = Theobald_import_exp[:, 1]
exp_dia = Theobald_import_exp[:, 2]
exp_bar = Theobald_import_exp[:, 3]
exp_deg = Theobald_import_exp[:, 4]

# define lists for iteration
max_range_out = []
max_height_out = []
max_height_dist_out = []

# iterate through data and build Plotting table
for runno in range(s_Theobald):
    linefile = chid_list[runno] + '_line.csv'
    out_import = pd.read_csv(outdir + linefile, skiprows=1)

    # import fds data
    AMPUAx = out_import['AMPUA-x'].values
    AMPUA = out_import['AMPUA'].values
    AMPUA = AMPUA[~np.isnan(AMPUA)]
    ZMAXz = out_import['ZMAX-z'].values
    ZMAX = out_import['ZMAX'].values
    ZMAX = ZMAX[~np.isnan(ZMAX)]
    ZMAX_Xz = out_import['ZMAX_X-z'].values
    ZMAX_X = out_import['ZMAX_X'].values

    # find FDS max range
    max_range = AMPUAx[-1]
    for f in range(len(AMPUA) - 1, -1, -1):
        max_range = AMPUAx[f]
        if AMPUA[f] > 1:
            break

    # find FDS max height
    max_height = ZMAXz[-1]
    for f in range(len(ZMAX) - 1, -1, -1):
        max_height = ZMAXz[f]
        if ZMAX[f] > 1:
            break

    # get max height dist
    max_height_dist = ZMAX_X[f]

    # build fds results lists
    max_range_out.append(max_range)
    max_height_out.append(max_height)
    max_height_dist_out.append(max_height_dist)

# Convert lists to numpy arrays
max_range_out = np.array(max_range_out)
max_height_out = np.array(max_height_out)
max_height_dist_out = np.array(max_height_dist_out)

# Scale Diameter for plots
exp_dia1 = exp_dia * 12

# Remove NAN rows from max height and max dist arrays, make tables
TheoArr1 = np.column_stack((max_range_exp, max_range_out, exp_noz, exp_dia1, exp_bar, exp_deg))
TheoArr2 = np.column_stack((max_height_exp, max_height_out, exp_noz, exp_dia1, exp_bar, exp_deg))
TheoArr3 = np.column_stack((max_height_dist_exp, max_height_dist_out, exp_noz, exp_dia1, exp_bar, exp_deg))
TheoArr2 = TheoArr2[~np.isnan(TheoArr2).any(axis=1)]
TheoArr3 = TheoArr3[~np.isnan(TheoArr3).any(axis=1)]
T1 = pd.DataFrame(TheoArr1, columns=['max range exp', 'max range out', 'nozzle', 'diameter', 'operating pressure', 'firing angle'])
T2 = pd.DataFrame(TheoArr2, columns=['max height exp', 'max height out', 'nozzle', 'diameter', 'operating pressure', 'firing angle'])
T3 = pd.DataFrame(TheoArr3, columns=['max height dist exp', 'max height dist out', 'nozzle', 'diameter', 'operating pressure', 'firing angle'])

# ----------------- Plot 1 Max Range -----------------

fig = fdsplotlib.plot_to_fig(x_data=[0,80], y_data=[0,80], marker_style='k-',
                             x_min=0, x_max=80, y_min=0, y_max=80,
                             figure_size=(plot_style['Scat_Paper_Width'],plot_style['Scat_Paper_Height']),
                             plot_size=(plot_style['Scat_Plot_Width'],plot_style['Scat_Plot_Height']),
                             plot_origin=(plot_style['Scat_Plot_X'],plot_style['Scat_Plot_Y']),
                             revision_label=version_string,
                             x_label='Measured Max Range (m)',
                             y_label='Predicted Max Range (m)')

range_6  = T1['nozzle'] == 6
range_7  = T1['nozzle'] == 7
range_9  = T1['nozzle'] == 9
range_10 = T1['nozzle'] == 10
fdsplotlib.plot_to_fig(x_data=T1.loc[range_7, 'max range exp'], y_data=T1.loc[range_7, 'max range out'], figure_handle=fig, marker_style='co', data_label='Nozzle 7')
fdsplotlib.plot_to_fig(x_data=T1.loc[range_9, 'max range exp'], y_data=T1.loc[range_9, 'max range out'], figure_handle=fig, marker_style='ko', data_label='Nozzle 9')
fdsplotlib.plot_to_fig(x_data=T1.loc[range_10,'max range exp'], y_data=T1.loc[range_10,'max range out'], figure_handle=fig, marker_style='mo', data_label='Nozzle 10')
fdsplotlib.plot_to_fig(x_data=T1.loc[range_6, 'max range exp'], y_data=T1.loc[range_6, 'max range out'], figure_handle=fig, marker_style='ro', data_label='Nozzle 6')

plt.savefig(plot_dir + 'Theobald_Hose_Stream_Max_Range.pdf', format='pdf')
plt.close(fig)

# ----------------- Plot 2 Max Height -----------------

fig = fdsplotlib.plot_to_fig(x_data=[0,20], y_data=[0,20], marker_style='k-',
                             x_min=0, x_max=20, y_min=0, y_max=20,
                             figure_size=(plot_style['Scat_Paper_Width'],plot_style['Scat_Paper_Height']),
                             plot_size=(plot_style['Scat_Plot_Width'],plot_style['Scat_Plot_Height']),
                             plot_origin=(plot_style['Scat_Plot_X'],plot_style['Scat_Plot_Y']),
                             revision_label=version_string,
                             x_label='Measured Max Height (m)',
                             y_label='Predicted Max Height (m)')

range_7 = T2['nozzle'] == 7
range_9 = T2['nozzle'] == 9
fdsplotlib.plot_to_fig(x_data=T2.loc[range_7, 'max height exp'], y_data=T2.loc[range_7, 'max height out'], figure_handle=fig, marker_style='co', data_label='Nozzle 7')
fdsplotlib.plot_to_fig(x_data=T2.loc[range_9, 'max height exp'], y_data=T2.loc[range_9, 'max height out'], figure_handle=fig, marker_style='ko', data_label='Nozzle 9')

plt.savefig(plot_dir + 'Theobald_Hose_Stream_Max_Height.pdf', format='pdf')
plt.close(fig)

# ----------------- Plot 3 Max Height Dist -----------------

fig = fdsplotlib.plot_to_fig(x_data=[0,50], y_data=[0,50], marker_style='k-',
                             x_min=0, x_max=50, y_min=0, y_max=50,
                             figure_size=(plot_style['Scat_Paper_Width'],plot_style['Scat_Paper_Height']),
                             plot_size=(plot_style['Scat_Plot_Width'],plot_style['Scat_Plot_Height']),
                             plot_origin=(plot_style['Scat_Plot_X'],plot_style['Scat_Plot_Y']),
                             revision_label=version_string,
                             x_label='Measured Max Height Distance (m)',
                             y_label='Predicted Max Height Distance (m)')

range_7 = T3['nozzle'] == 7
range_9 = T3['nozzle'] == 9
fdsplotlib.plot_to_fig(x_data=T3.loc[range_7, 'max height dist exp'], y_data=T3.loc[range_7, 'max height dist out'], figure_handle=fig, marker_style='co', data_label='Nozzle 7')
fdsplotlib.plot_to_fig(x_data=T3.loc[range_9, 'max height dist exp'], y_data=T3.loc[range_9, 'max height dist out'], figure_handle=fig, marker_style='ko', data_label='Nozzle 9')

plt.savefig(plot_dir + 'Theobald_Hose_Stream_Max_Height_Distance.pdf', format='pdf')
plt.close(fig)

