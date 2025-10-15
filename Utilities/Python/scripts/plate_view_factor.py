
# Plate View Factor

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import fdsplotlib

plot_style = fdsplotlib.get_plot_style('fds')

outdir = '../../Verification/Radiation/'
pltdir = '../../Manuals/FDS_Verification_Guide/SCRIPT_FIGURES/'

git_file = outdir + 'plate_view_factor_2D_30_git.txt'
version_string = fdsplotlib.get_version_string(git_file)

NRA = np.array([30, 60, 100])

Exact_Flux_2D = 105.34
Exact_Flux_cart = 81.8
Exact_Flux_cyl = 74.1

Flux_2D = np.zeros(3)
Flux_cart = np.zeros(3)
Flux_cyl = np.zeros(3)
Flux_ibm = np.zeros(3)

# 2D
M = pd.read_csv(outdir+'plate_view_factor_2D_30_devc.csv', skiprows=2, header=None)
Flux_2D[0] = M.iloc[:, 1].max()
M = pd.read_csv(outdir+'plate_view_factor_2D_60_devc.csv', skiprows=2, header=None)
Flux_2D[1] = M.iloc[:, 1].max()
M = pd.read_csv(outdir+'plate_view_factor_2D_100_devc.csv', skiprows=2, header=None)
Flux_2D[2] = M.iloc[:, 1].max()

# Cart
M = pd.read_csv(outdir+'plate_view_factor_cart_30_devc.csv', skiprows=2, header=None)
Flux_cart[0] = M.iloc[:, 1].max()
M = pd.read_csv(outdir+'plate_view_factor_cart_60_devc.csv', skiprows=2, header=None)
Flux_cart[1] = M.iloc[:, 1].max()
M = pd.read_csv(outdir+'plate_view_factor_cart_100_devc.csv', skiprows=2, header=None)
Flux_cart[2] = M.iloc[:, 1].max()

# Cylindrical
M = pd.read_csv(outdir+'plate_view_factor_cyl_30_devc.csv', skiprows=2, header=None)
Flux_cyl[0] = M.iloc[:, 1].max()
M = pd.read_csv(outdir+'plate_view_factor_cyl_60_devc.csv', skiprows=2, header=None)
Flux_cyl[1] = M.iloc[:, 1].max()
M = pd.read_csv(outdir+'plate_view_factor_cyl_100_devc.csv', skiprows=2, header=None)
Flux_cyl[2] = M.iloc[:, 1].max()

# IBM
M = pd.read_csv(outdir+'plate_view_factor_ibm_30_devc.csv', skiprows=2, header=None)
Flux_ibm[0] = M.iloc[:, 1].max()
M = pd.read_csv(outdir+'plate_view_factor_ibm_60_devc.csv', skiprows=2, header=None)
Flux_ibm[1] = M.iloc[:, 1].max()
M = pd.read_csv(outdir+'plate_view_factor_ibm_100_devc.csv', skiprows=2, header=None)
Flux_ibm[2] = M.iloc[:, 1].max()


fig = fdsplotlib.plot_to_fig(x_data=[0,0], y_data=[0,0],
                             x_min=20, x_max=110, y_min=30, y_max=120,
                             revision_label=version_string,
                             plot_title=r'Radiative Heat Flux (plate_view_factor)',
                             x_label='Number of Radiation Angles',
                             y_label='Heat Flux (kW/m$^2$)')

fdsplotlib.plot_to_fig(x_data=NRA, y_data=Exact_Flux_2D*np.ones(3), marker_style='r-', figure_handle=fig)
fdsplotlib.plot_to_fig(x_data=NRA, y_data=Flux_2D                 , marker_style='ro', data_label='FDS 2D', figure_handle=fig)
fdsplotlib.plot_to_fig(x_data=NRA, y_data=Exact_Flux_cart*np.ones(3), marker_style='b-', figure_handle=fig)
fdsplotlib.plot_to_fig(x_data=NRA, y_data=Flux_cart                 , marker_style='bs', data_label='FDS 3D', figure_handle=fig)
fdsplotlib.plot_to_fig(x_data=NRA, y_data=Exact_Flux_cyl*np.ones(3), marker_style='g-', figure_handle=fig)
fdsplotlib.plot_to_fig(x_data=NRA, y_data=Flux_cyl                 , marker_style='gd', data_label='FDS Cyl.', figure_handle=fig)

plt.savefig(pltdir + 'plate_view_factor.pdf', format='pdf')
plt.close(fig)

