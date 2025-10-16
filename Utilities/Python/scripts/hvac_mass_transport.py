#!/usr/bin/env python3
# B M Ralph
# 12-8-2016
# hvac_mass_transport.py
#
# Converted by Floyd
# 10-14-2025
#
# Convergence study for HVAC transient mass transport (mass fraction at
# downstream duct node).

import os
import numpy as np
import pandas as pd
import math
import matplotlib.pyplot as plt
import fdsplotlib


import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

datadir =  '../../Verification/HVAC/'
plotdir = '../../Manuals/FDS_Verification_Guide/SCRIPT_FIGURES/';
filename = [
    'HVAC_mass_transport_conv_0020_devc.csv',
    'HVAC_mass_transport_conv_0040_devc.csv',
    'HVAC_mass_transport_conv_0080_devc.csv',
    'HVAC_mass_transport_conv_0160_devc.csv',
    'HVAC_mass_transport_conv_0320_devc.csv'
]
PlotStyle = ['b-', 'g-', 'r-', 'm-', 'c-']
nc_array = np.array([20, 40, 80, 160, 320])
dx_array = np.array([1/20, 1/40, 1/80, 1/160, 1/320])

# Input parameters
t0 = 0.0
t_end = 2.0
u = 1.0
L = 1.0

skip_case = False

for i in range(len(filename)):
   name = datadir+filename[i]
   if not os.path.exists(name):
      skip_case = True
      print('Error: File ', filename[i], ' does not exist. Skipping case.')

if skip_case: quit()

# --- Analytical Solution ---

# Analytical solution (using a lambda function which is similar to MATLAB's Anonymous Function)
# Y = t * (t <= 1) + 1 * (t > 1)
# Note: numpy arrays are required for element-wise comparison and multiplication

# Create mass fraction plot and plot analytical solution
nt = 1000
dt = t_end / nt
# Generate time vector similar to MATLAB's t0 : dt : t_end
tc = np.arange(t0, t_end + dt, dt)
Y =(lambda t: (t * (t <= 1) + 1.0 * (t > 1)).astype(float))(tc)

L1e = [] # initialize L1 norm error vector
L2e = [] # initialize L2 norm error vector
dxx = [] # init dxx vector

git_file = datadir+'HVAC_mass_transport_conv_0320_git.txt'
version_string = fdsplotlib.get_version_string(git_file)

fig = fdsplotlib.plot_to_fig(x_data=tc, y_data=Y, marker_style='k--',
      revision_label=version_string,x_min=0,x_max=2,y_min=0,y_max=1.2,
      data_label = 'Exact Solution',
      x_label='Time (s)',
      y_label='Mass Fraction (kg/kg)')

for i in range(len(filename)):
   df = pd.read_csv(datadir+filename[i], skiprows=2) # The header is in row 1 (0-indexed)

   # Get column names (assuming standard FDS output of Time, Variable 1, Variable 2, ...)
   # The MATLAB code uses M.data(1:end,1) for Time and M.data(1:end,2) for the species data.
   t_fds = df.iloc[:, 0].values # Time (1st column, 0-indexed)
   Y_fds = df.iloc[:, 1].values # FDS species data (2nd column, 1-indexed)

   # 3. Plot FDS results
   fdsplotlib.plot_to_fig(x_data=t_fds, y_data=Y_fds, marker_style=PlotStyle[i],
   figure_handle=fig,
   data_label='FDS N_CELLS = '+str(nc_array[i]))

   # 4. Compute Error Norms (element-by-element appending)

   # Evaluate Analytical solution at FDS time points
   Y_analytical_at_t_fds = (lambda t: (t * (t <= 1) + 1.0 * (t > 1)).astype(float))(t_fds)

   # L1 norm error: 1/nc * sum(|Y_analytical - Y_fds|)
   L1_error = (1 / nc_array[i]) * np.sum(np.abs(Y_analytical_at_t_fds - Y_fds))
   L1e.append(L1_error)

   # L2 norm error: sqrt(1/nc * sum((Y_analytical - Y_fds)^2))
   L2_error = np.sqrt((1 / nc_array[i]) * np.sum((Y_analytical_at_t_fds - Y_fds)**2))
   L2e.append(L2_error)

plotname = plotdir + 'HVAC_mass_transport_convergence_1.pdf'
plt.savefig(plotname, format='pdf')
plt.close()

fig = fdsplotlib.plot_to_fig(x_data=dx_array, y_data=L1e, marker_style='ks-',plot_type='loglog',
      revision_label=version_string,x_min=0.001,x_max=0.1,y_min=1E-4,y_max=0.1,
      data_label = 'FDS',
      x_label='$\\Delta x$ (m)',
      y_label='L$_1$ error')

fdsplotlib.plot_to_fig(x_data=dx_array, y_data=dx_array/10, marker_style='k--',
figure_handle=fig,
data_label='$O$($\\Delta x$)')


plotname = plotdir + 'HVAC_mass_transport_convergence_2.pdf'
plt.savefig(plotname, format='pdf')
plt.close()

fig = fdsplotlib.plot_to_fig(x_data=dx_array, y_data=L2e, marker_style='ks-',plot_type='loglog',
      revision_label=version_string,x_min=0.001,x_max=0.1,y_min=1E-4,y_max=0.1,
      data_label = 'FDS',
      x_label='$\\Delta x$ (m)',
      y_label='L$_2$ error')

fdsplotlib.plot_to_fig(x_data=dx_array, y_data=dx_array/10, marker_style='k--',
figure_handle=fig,
data_label='$O$($\\Delta x$)')


plotname = plotdir + 'HVAC_mass_transport_convergence_3.pdf'
plt.savefig(plotname, format='pdf')
plt.close()

