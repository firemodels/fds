
# Plot Werner and Wengle velocity profile for FDS Technical Reference Guide

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import os
import fdsplotlib

plot_style = fdsplotlib.get_plot_style('fds')

pltdir = '../../Manuals/FDS_Technical_Reference_Guide/SCRIPT_FIGURES/'

A = 8.3
B = 1/7
n = 100
zp = np.logspace(0, 4, n)
uu = zp.copy()  # log law
up = np.zeros(n)

for j in range(n):
    if zp[j] < 11.81:
        up[j] = zp[j]
    else:
        up[j] = A * zp[j]**B  # Werner and Wengle (power law)
        uu[j] = (1/0.41) * np.log(zp[j]) + 5.2  # log law

# Create figure
fig = fdsplotlib.plot_to_fig(x_data=[11.8,11.8], y_data=[0,32], marker_style='k--',
                             x_min=1, x_max=10000, y_min=0, y_max=32,
                             plot_type='semilogx',
                             x_label='$z^+$',
                             y_label='$u^+$')

fdsplotlib.plot_to_fig(x_data=zp, y_data=up, figure_handle=fig, marker_style='b-')
fdsplotlib.plot_to_fig(x_data=zp, y_data=uu, figure_handle=fig, marker_style='r--')

# Add text annotations
ax = plt.gca()
ax.text(1.5, 8, r'$u^+ = z^+$', fontsize=plot_style['Label_Font_Size'])
ax.text(200, 15, r'$u^+ = 2.4 \ln z^+ + 5.2$', fontsize=plot_style['Label_Font_Size'], color=[1, 0, 0])
ax.text(500, 30, r'$u^+ = A(z^+)^B$', fontsize=plot_style['Label_Font_Size'], color=[0, 0, 1])
ax.text(15, 5, r'$z^+ = 11.81$', fontsize=plot_style['Label_Font_Size'])

fig.savefig(pltdir + 'lawofthewall.pdf', format='pdf')
plt.close()

