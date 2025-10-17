# McGrattan
# 2-16-15
# radiating_polygon.m
#
# Converted by Floyd
# 10/16/2025

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import fdsplotlib

datadir = '../../Verification/Radiation/'
plotdir = '../../Manuals/FDS_Verification_Guide/SCRIPT_FIGURES/'
filename = ['radiating_polygon_square_20_line.csv','radiating_polygon_square_40_line.csv','radiating_polygon_square_80_line.csv']

skip_case = False

for i in range(len(filename)):
   name = datadir+filename[i]
   if not os.path.exists(name):
      skip_case = True
      print('Error: File ', filename[i], ' does not exist. Skipping case.')

if skip_case: quit()

flux_20 = pd.read_csv(datadir+filename[0],skiprows=2,header=None)
flux_40 = pd.read_csv(datadir+filename[1],skiprows=2,header=None)
flux_80 = pd.read_csv(datadir+filename[2],skiprows=2,header=None)

n = 4
pi = np.pi

z_max = 1
z_min = 1 - 0.001 * (1000 - 1)  # z(1000) = 1 - 0.001 * 999 = 0.001

z = np.linspace(z_max, z_min, 1000)

R = 0.5 * np.sqrt(2) / z
H = 0.5 / z

# flux(j)=5.67e-11*1273.15^4*n*H/(pi*sqrt(1+H^2))*atan(sqrt((R^2-H^2)/(1+H^2)));
T_hot = 1273.15
sigma = 5.67e-11

term1 = sigma * (T_hot**4) * n
term2 = H / (pi * np.sqrt(1 + H**2))
term3 = np.arctan(np.sqrt((R**2 - H**2) / (1 + H**2)))

flux = term1 * term2 * term3

git_file = datadir+'radiating_polygon_square_20_git.txt'
version_string = fdsplotlib.get_version_string(git_file)

fig = fdsplotlib.plot_to_fig(x_data=z, y_data=flux, marker_style='k-',
      revision_label=version_string,x_min=0,x_max=1,y_min=20,y_max=160,
      plot_title = 'Radiative Flux from a Hot Square Plate',
      data_label='Exact',
      y_label='Radiative Heat Flux (kW/m$^2$)',
      x_label='Distance from Plate (m)')

fdsplotlib.plot_to_fig(x_data=1-flux_20[0], y_data=flux_20[1], marker_style='r-',
      figure_handle=fig,
      data_label='5 cm, 100 angles')

fdsplotlib.plot_to_fig(x_data=1-flux_40[0], y_data=flux_40[1], marker_style='g-',
      figure_handle=fig,
      data_label='2.5 cm, 400 angles')

fdsplotlib.plot_to_fig(x_data=1-flux_80[0], y_data=flux_80[1], marker_style='b-',
      figure_handle=fig,
      data_label='1.25 cm, 1000 angles')

plotname = plotdir + 'radiating_polygon_square.pdf'
plt.savefig(plotname, format='pdf')
plt.close()
