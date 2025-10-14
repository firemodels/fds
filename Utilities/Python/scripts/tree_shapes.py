import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import fdsplotlib

plot_style = fdsplotlib.get_plot_style('fds')

# Paths and identifiers
outdir = '../../Verification/WUI/'
pltdir = '../../Manuals/FDS_Verification_Guide/SCRIPT_FIGURES/'
CHID = 'tree_shapes'

devc_file = os.path.join(outdir, f'{CHID}_devc.csv')
if not os.path.exists(devc_file):
    print(f'Error: File {devc_file} does not exist. Skipping case.')
    exit(1)

# Read device data
DEVC = pd.read_csv(devc_file, skiprows=1)

# Extract data from final row
random_part_mass = DEVC.iloc[-1][1:6].values
one_part_mass =  DEVC.iloc[-1][6:].values

# Ideal input volumes
box_vol = 0.6 * 0.6 * 0.4
cone_vol = 0.6 * (0.3 ** 2) * np.pi / 3.0
cyl_vol = 0.6 * (0.2 ** 2) * np.pi
cone_shell_vol = 0.6 * ((0.3 ** 2) - (0.2 ** 2)) * np.pi / 3.0
cyl_shell_vol = 0.6 * ((0.2 ** 2) - (0.15 ** 2)) * np.pi

input_vol = np.array([box_vol, cone_vol, cyl_vol, cone_shell_vol, cyl_shell_vol])

# Relative error in mass (bulk density 5 kg/m^3)
expected_mass = 5.0 * np.concatenate([input_vol, input_vol])
observed_mass = np.concatenate([random_part_mass, one_part_mass])
rel_err = np.abs(observed_mass - expected_mass) / expected_mass
rel_err = np.nan_to_num(rel_err,nan=1.)
if np.max(rel_err) > 0.1:
    print('Python Warning: The mass in tree_shapes is out of tolerance.')

# Revision string
git_file = os.path.join(outdir, f'{CHID}_git.txt')
version_string = fdsplotlib.get_version_string(git_file)

# Plot expected value
fig = fdsplotlib.plot_to_fig(
    x_data=[0.0, 0.2],
    y_data=[0.0, 1.0],
    marker_style='k--',
    x_min=0.0,
    x_max=0.15,
    y_min=0.0,
    y_max=0.8,
    data_label='ideal',
    revision_label=version_string,
    x_label=r'Input volume (m$^3$)',
    y_label='Tree crown mass (kg)'
)

# Plot results
fdsplotlib.plot_to_fig(
    x_data=input_vol,
    y_data=random_part_mass,
    marker_style='ro',
    data_label='1000 random particles',
    figure_handle=fig)
fdsplotlib.plot_to_fig(
    x_data=input_vol,
    y_data=one_part_mass,
    marker_style='bo',
    data_label='1 particle per cell',
    figure_handle=fig)

# Save figure
fig.savefig(os.path.join(pltdir, 'tree_shapes.pdf'), format='pdf')
plt.close(fig)

