
# Makes figure for HVAC Fan Parameters section of the FDS Users Guide

import numpy as np
import matplotlib.pyplot as plt
import fdsplotlib

plot_style = fdsplotlib.get_plot_style('fds')

pltdir = '../../Manuals/FDS_User_Guide/SCRIPT_FIGURES/'

vdot_max = 10
dp_max = 500
dp = np.arange(-1000, 1001, 1)

# Constant volume flow rate
vdot1 = 10

# Plot quadratic fan curve (black line)
vdot = vdot_max * np.sign(dp_max - dp) * np.sqrt(np.abs(dp - dp_max) / dp_max)

# Create user fan curve (ramp) data points
rampx = []
rampy = []
i = 0
for dp_val in range(-1000, 1001, 200):
    rampx.append(vdot_max * np.sign(dp_max - dp_val) * np.sqrt(np.abs(dp_val - dp_max) / dp_max))
    rampy.append(dp_val)
    i = i + 1

fig = fdsplotlib.plot_to_fig(x_data=[vdot1,vdot1], y_data=[-1000,1000], marker_style='r-', data_label='constant volume',
                             x_min=-10, x_max=20, y_min=-1000, y_max=1000,
                             x_label='Volume Flow Rate (m$^3$/s)',
                             y_label='Static Pressure (Pa)')

fdsplotlib.plot_to_fig(x_data=vdot, y_data=dp, marker_style='k-', data_label='quadratic', figure_handle=fig)
fdsplotlib.plot_to_fig(x_data=rampx, y_data=rampy, marker_style='b-', data_label='user fan curve', figure_handle=fig)

# Enable grid
ax = plt.gca()
ax.grid(True, which='both', axis='both')

plt.savefig(pltdir + 'fan_curve.pdf', format='pdf')
plt.close()

