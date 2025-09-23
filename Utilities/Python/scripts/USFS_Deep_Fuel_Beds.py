
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.integrate import trapezoid

# include FDS plot styles, etc.
import fdsplotlib

# Get plot style parameters
plot_style = fdsplotlib.get_plot_style('fds')

# Paths
outdir = '../../../out/USFS_Deep_Fuel_Beds/'
expdir = '../../../exp/USFS_Deep_Fuel_Beds/'
plot_dir = '../../Manuals/FDS_Validation_Guide/SCRIPT_FIGURES/USFS_Deep_Fuel_Beds/'

git_file =  os.path.join(outdir, f'burn100_60D_15S_35L_cat_git.txt')
version_string = fdsplotlib.get_version_string(git_file)
    
# Read experimental data
EXP = pd.read_csv(os.path.join(expdir, 'exp_params.csv'), header=0)
BURN_NO = EXP['BURN_NO'].values
SPACING = EXP['SPACING'].values * 2.5
DEPTH = EXP['DEPTH'].values * 2.5
FMC = EXP['MOISTURE_CONTENT'].values
SLOPE = EXP['SLOPE'].values
BURN_PCT_EXP = EXP['BURN_PERCENT'].values
BURN_PCT_FDS = np.zeros(len(BURN_NO))
    
for i in range(len(BURN_NO)):
    CHID = f'burn{int(BURN_NO[i])}_{int(DEPTH[i])}D_{int(SLOPE[i])}S_{int(SPACING[i])}L'
    HRR = pd.read_csv(os.path.join(outdir, f'{CHID}_cat_hrr.csv'), header=1)
    HRR.columns = HRR.columns.str.strip()
    DEVC = pd.read_csv(os.path.join(outdir, f'{CHID}_cat_devc.csv'), header=1)
    DEVC.columns = DEVC.columns.str.strip()
    # Determine FDS (out) Burn Type:
    hrr_col = HRR['HRR'].values  # (kW)
    time_col = HRR['Time'].values  # (s)
    # Note: adjust column name if actual CSV has quotes around "Fuel Load"
    fuel_load = DEVC.iloc[1].get('Fuel Load', DEVC.iloc[1].get('"Fuel Load"', np.nan))
    dry_fuel_mass = (1 / (1 + FMC[i] * 0.01)) * fuel_load  # (kg)
    total_fuel_energy = dry_fuel_mass * 17425.0  # (kJ)
    ignition_energy = trapezoid([500 * 0.3 * 1.6, 500 * 0.3 * 1.6, 0], [0, 10, 20])  # (kJ)
    total_fds_energy = trapezoid(hrr_col, time_col)  # (kJ)
    BURN_PCT_FDS[i] = (total_fds_energy - ignition_energy) / total_fuel_energy * 100.0

fig = fdsplotlib.plot_to_fig([0,100],[0,100],
                             x_min=0, x_max=100, y_min=0, y_max=100,
                             figure_size=(plot_style['Scat_Paper_Width'],plot_style['Scat_Paper_Height']),
                             plot_size=(plot_style['Scat_Plot_Width'],plot_style['Scat_Plot_Height']),
                             plot_origin=(plot_style['Scat_Plot_X'],plot_style['Scat_Plot_Y']),
                             revision_label=version_string,
                             x_label='Measured Fuel Consumption (\%)',
                             y_label='Predicted Fuel Consumption (\%)'
                             )

fdsplotlib.plot_to_fig(BURN_PCT_EXP, BURN_PCT_FDS, figure_handle=fig, marker_style='ko')
    
for i in range(len(BURN_NO)):
    plt.text(BURN_PCT_EXP[i], BURN_PCT_FDS[i], f'{"  "}+{int(BURN_NO[i])}', fontname=plot_style['Font_Name'], fontsize=6)

plt.savefig(os.path.join(plot_dir, 'USFS_Deep_Fuel_Beds_Scatterplot.pdf'), format='pdf')


