
# Postprocessing script to extract transformation mapping from SMV and create figures for User's Guide

import numpy as np
import matplotlib.pyplot as plt
import fdsplotlib

plot_style = fdsplotlib.get_plot_style('fds')

outdir = '../../Verification/Miscellaneous/'
pltdir = '../../Manuals/FDS_User_Guide/SCRIPT_FIGURES/'

fid = open(outdir + 'mesh_transformation.smv', 'r')

I = 50
J = 50
LX = 1.5
LY = 1.5
dx = LX / I
dy = LY / J

cx = np.zeros(I + 1)
x = np.zeros(I + 1)
cy = np.zeros(J + 1)
y = np.zeros(J + 1)

while True:
    tline = fid.readline()
    if not tline:  # End of file
        break
    tline = tline.strip()
    if tline == 'TRNX':
        # Skip next 3 lines
        for i in range(3):
            skip = fid.readline()
        # Read I+1 lines of transformation data
        for i in range(I + 1):
            tline = fid.readline()
            values = [float(val) for val in tline.split()]
            cx[i] = values[0] * dx
            x[i] = values[1]
    if tline == 'TRNY':
        # Skip next 3 lines
        for j in range(3):
            skip = fid.readline()
        # Read J+1 lines of transformation data
        for j in range(J + 1):
            tline = fid.readline()
            values = [float(val) for val in tline.split()]
            cy[j] = values[0] * dy
            y[j] = values[1]

fid.close()

# From the mesh_transformation.fds input file:
CC = [0.3, 1.2]
PC = [0.5, 1.0]

XMinorTick = np.arange(0.1, 1.5, 0.1)

# First figure - piece-wise linear transformation

fig = fdsplotlib.plot_to_fig(x_data=cx, y_data=x, marker_style='k-',
                             x_min=0, x_max=1.5, y_min=0, y_max=1.5,
                             figure_size=(plot_style['Scat_Paper_Width'],plot_style['Scat_Paper_Height']),
                             plot_size=(plot_style['Scat_Plot_Width'],plot_style['Scat_Plot_Height']),
                             plot_origin=(plot_style['Scat_Plot_X'],plot_style['Scat_Plot_Y']),
                             x_label=r'$\xi$', y_label=r'$x$')

for i in range(len(XMinorTick)):
    if XMinorTick[i] < CC[0]:
        YMinorTick = (PC[0] / CC[0]) * XMinorTick[i]
    elif XMinorTick[i] < CC[1]:
        YMinorTick = PC[0] + (PC[1] - PC[0]) / (CC[1] - CC[0]) * (XMinorTick[i] - CC[0])
    elif XMinorTick[i] <= LX:
        YMinorTick = PC[1] + (LX - PC[1]) / (LX - CC[1]) * (XMinorTick[i] - CC[1])

    fdsplotlib.plot_to_fig(x_data=[XMinorTick[i],XMinorTick[i]], y_data=[0,YMinorTick], marker_style='k-', figure_handle=fig)
    fdsplotlib.plot_to_fig(x_data=[0, XMinorTick[i]], y_data=[YMinorTick,YMinorTick],   marker_style='k-', figure_handle=fig)

ax = plt.gca()
ax.set_xticks(np.arange(0, 1.5, 0.3))
ax.set_yticks(np.arange(0, 1.5, 0.3))

plt.savefig(pltdir + 'piece_wise_linear_trnx.pdf', format='pdf')
plt.close()

# Second figure - polynomial transformation

fig = fdsplotlib.plot_to_fig(x_data=cy, y_data=y, marker_style='k-',
                             x_min=0, x_max=1.5, y_min=0, y_max=1.5,
                             figure_size=(plot_style['Scat_Paper_Width'],plot_style['Scat_Paper_Height']),
                             plot_size=(plot_style['Scat_Plot_Width'],plot_style['Scat_Plot_Height']),
                             plot_origin=(plot_style['Scat_Plot_X'],plot_style['Scat_Plot_Y']),
                             x_label=r'$\xi$', y_label=r'$x$')

for i in range(len(XMinorTick)):
    YMinorTick = 2 * XMinorTick[i] - 2 * XMinorTick[i]**2 + 0.8889 * XMinorTick[i]**3
    fdsplotlib.plot_to_fig(x_data=[XMinorTick[i],XMinorTick[i]], y_data=[0,YMinorTick], marker_style='k-', figure_handle=fig)
    fdsplotlib.plot_to_fig(x_data=[0, XMinorTick[i]], y_data=[YMinorTick,YMinorTick],   marker_style='k-', figure_handle=fig)

ax = plt.gca()
ax.set_xticks(np.arange(0, 1.5, 0.3))
ax.set_yticks(np.arange(0, 1.5, 0.3))

plt.savefig(pltdir + 'polynomial_trnx.pdf', format='pdf')
plt.close()

