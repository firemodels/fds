
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
import os
import fdsplotlib

plot_style = fdsplotlib.get_plot_style('fds')

outdir = '../../../out/NIST_NRC_Parallel_Panels/'
expdir = '../../../exp/Submodules/macfp-db/Fire_Growth/NIST_Parallel_Panel/Experimental_Data/'
pltdir = '../../Manuals/FDS_Validation_Guide/SCRIPT_FIGURES/NIST_NRC_Parallel_Panels/'

git_file = outdir + 'PMMA_60_kW_1_cm_git.txt'
version_string = fdsplotlib.get_version_string(git_file)

color = ['r','g','c','k','m','b','y']
HRR_label = ['120 kW', '200 kW', '300 kW', '500 kW', '1500 kW', '2000 kW', '2800 kW']
Time_label = ['20 s', '40 s', '60 s', '80 s']

def importdata(filename, delimiter=',', headerlines=0):
    """
    Mimics MATLAB's importdata function for CSV files
    Returns a structure with .data attribute containing numeric data
    """
    class DataStruct:
        def __init__(self):
            self.data = None

    result = DataStruct()
    result.data = pd.read_csv(filename, delimiter=delimiter, skiprows=headerlines, header=None).values
    return result

# First Plot: Vertical Heat Flux for PMMA Fire

z = np.array([10, 20, 30, 50, 75, 100, 140, 180, 220])
zm = np.linspace(1, 243, 50)

E  = importdata(expdir + 'PMMA_heatflux.csv', ',', 2)
M4 = importdata(outdir + 'PMMA_60_kW_4_cm_devc.csv', ',', 2)
M2 = importdata(outdir + 'PMMA_60_kW_2_cm_devc.csv', ',', 2)
M1 = importdata(outdir + 'PMMA_60_kW_1_cm_devc.csv', ',', 2)
H4 = importdata(outdir + 'PMMA_60_kW_4_cm_hrr.csv', ',', 2)
H2 = importdata(outdir + 'PMMA_60_kW_2_cm_hrr.csv', ',', 2)
H1 = importdata(outdir + 'PMMA_60_kW_1_cm_hrr.csv', ',', 2)

indices4 = []
indices2 = []
indices1 = []
j = 0
for i in [1, 2, 3, 5, 8, 9, 10]:
    # For M4
    MM = np.max(H4.data[:, 1])
    II = np.argmax(H4.data[:, 1])
    found_indices = np.where(H4.data[:, 1] > E.data[i-1, 0])[0]
    if len(found_indices) > 0:
        indices4.append(min(np.min(found_indices), II))
    else:
        indices4.append(II)

    # For M2
    MM = np.max(H2.data[:, 1])
    II = np.argmax(H2.data[:, 1])
    found_indices = np.where(H2.data[:, 1] > E.data[i-1, 0])[0]
    if len(found_indices) > 0:
        indices2.append(min(np.min(found_indices), II))
    else:
        indices2.append(II)

    # For M1
    MM = np.max(H1.data[:, 1])
    II = np.argmax(H1.data[:, 1])
    found_indices = np.where(H1.data[:, 1] > E.data[i-1, 0])[0]
    if len(found_indices) > 0:
        indices1.append(min(np.min(found_indices), II))
    else:
        indices1.append(II)
    j = j + 1

fig = fdsplotlib.plot_to_fig(x_data=[-1,-1], y_data=[-1,-1],
                             x_min=0, x_max=150, y_min=0, y_max=250,
                             revision_label=version_string,
                             x_label='Heat Flux (kW/m$^2$)',
                             y_label='Height (cm)')

qdot = {}
err = {}
qdotm = {}

# Plot experimental data with error bars
color_idx = 0
for i in [1, 2, 3, 5, 8, 9, 10]:
    qdot[i] = E.data[i-1, 1:10]
    err[i] = E.data[i-1, 10:19]
    fdsplotlib.plot_to_fig(x_data=qdot[i], y_data=z, x_error=err[i], figure_handle=fig, marker_style=color[color_idx]+'o', data_label=HRR_label[color_idx])
    color_idx += 1

# Plot M4 data
color_idx = 0
for idx, i in enumerate([1, 2, 3, 5, 8, 9, 10]):
    qdotm[indices4[idx]] = M4.data[indices4[idx], 1:51]
    fdsplotlib.plot_to_fig(x_data=qdotm[indices4[idx]], y_data=zm, figure_handle=fig, marker_style=color[color_idx]+'--')
    color_idx += 1

# Plot M2 data
color_idx = 0
for idx, i in enumerate([1, 2, 3, 5, 8, 9, 10]):
    qdotm[indices2[idx]] = M2.data[indices2[idx], 1:51]
    fdsplotlib.plot_to_fig(x_data=qdotm[indices2[idx]], y_data=zm, figure_handle=fig, marker_style=color[color_idx]+'-.')
    color_idx += 1

# Plot M1 data
color_idx = 0
for idx, i in enumerate([1, 2, 3, 5, 8, 9, 10]):
    qdotm[indices1[idx]] = M1.data[indices1[idx], 1:51]
    fdsplotlib.plot_to_fig(x_data=qdotm[indices1[idx]], y_data=zm, figure_handle=fig, marker_style=color[color_idx]+'-')
    color_idx += 1

fig.savefig(pltdir + 'PMMA_Heat_Flux.pdf', format='pdf')
plt.close(fig)


del z, E, M1, M2, qdot

# Second plot - Marinite Centerline Heat Flux

z = np.array([20, 40, 60, 80])
E  = importdata(expdir + 'Burner_HF_Centerline_multi-layer.csv', ',', 2)
M4 = importdata(outdir + 'Marinite_60_kW_4_cm_devc.csv', ',', 2)
M2 = importdata(outdir + 'Marinite_60_kW_2_cm_devc.csv', ',', 2)
M1 = importdata(outdir + 'Marinite_60_kW_1_cm_devc.csv', ',', 2)

fig = fdsplotlib.plot_to_fig(x_data=[-1,-1], y_data=[-1,-1],
                             x_min=0, x_max=80, y_min=0, y_max=180,
                             revision_label=version_string,
                             x_label='Heat Flux (kW/m$^2$)',
                             y_label='Height (cm)')

# Plot experimental data
color_idx = 0
for i in [1, 2, 3, 4]:
    qdot = {}
    qdot[i] = E.data[i-1, 2:6]
    fdsplotlib.plot_to_fig(x_data=qdot[i], y_data=z, figure_handle=fig, marker_style=color[color_idx]+'o', data_label=Time_label[color_idx])
    color_idx += 1

# Plot M4 data
color_idx = 0
for i in [15, 21, 31, 41]:
    fdsplotlib.plot_to_fig(x_data=M4.data[i-1, 1:51], y_data=zm, figure_handle=fig, marker_style=color[color_idx]+'-.')
    color_idx += 1

# Plot M2 data
color_idx = 0
for i in [15, 21, 31, 41]:
    fdsplotlib.plot_to_fig(x_data=M2.data[i-1, 1:51], y_data=zm, figure_handle=fig, marker_style=color[color_idx]+'--')
    color_idx += 1

# Plot M1 data
color_idx = 0
for i in [15, 21, 31, 41]:
    fdsplotlib.plot_to_fig(x_data=M1.data[i-1, 1:51], y_data=zm, figure_handle=fig, marker_style=color[color_idx]+'-')
    color_idx += 1

ax = plt.gca()
ax.set_yticks([0,30,60,90,120,150,180])
ax.text(78,110, 'Solid: 1 cm', horizontalalignment='right', fontsize=plot_style['Key_Font_Size'])
ax.text(78,100, 'Dashed: 2 cm', horizontalalignment='right', fontsize=plot_style['Key_Font_Size'])
ax.text(78, 90, 'Dash-Dot: 4 cm', horizontalalignment='right', fontsize=plot_style['Key_Font_Size'])

fig.savefig(pltdir + 'Marinite_Heat_Flux.pdf', format='pdf')
plt.close(fig)


# Third plot - Marinite Heat Flux Contours

del E
E = importdata(expdir + 'Burner_steadyHF_Width_multi-layer.csv', ',', 2)

x = np.array([-25, -15, 0, 15, 25])
y = np.array([20, 50, 75, 100])
X, Y = np.meshgrid(x, y)

Z = np.zeros((4, 5))
for i in range(4):
    Z[i, 0:5] = E.data[i, 1:6]

Z4 = np.zeros((4, 5))
Z4[0, 0:5] = M4.data[-1, 101:106]
Z4[1, 0:5] = M4.data[-1, 106:111]
Z4[2, 0:5] = M4.data[-1, 111:116]
Z4[3, 0:5] = M4.data[-1, 116:121]

Z2 = np.zeros((4, 5))
Z2[0, 0:5] = M2.data[-1, 101:106]
Z2[1, 0:5] = M2.data[-1, 106:111]
Z2[2, 0:5] = M2.data[-1, 111:116]
Z2[3, 0:5] = M2.data[-1, 116:121]

Z1 = np.zeros((4, 5))
Z1[0, 0:5] = M1.data[-1, 101:106]
Z1[1, 0:5] = M1.data[-1, 106:111]
Z1[2, 0:5] = M1.data[-1, 111:116]
Z1[3, 0:5] = M1.data[-1, 116:121]

levels = [5, 10, 15, 20, 30, 40, 50, 60]

fig = fdsplotlib.plot_to_fig(x_data=[-1,-1], y_data=[-1,-1],
                             x_min=-30, x_max=30, y_min=0, y_max=110,
                             revision_label=version_string,
                             x_label='Width (cm)',
                             y_label='Height (cm)')

# Filled contour for experimental data
C_exp = plt.contourf(X, Y, Z, levels=levels, linestyles='-', cmap='viridis')
plt.clabel(C_exp, fontsize=6, colors='black')

# Contour lines for FDS 4 cm
C4 = plt.contour(X, Y, Z4, levels=levels, colors='red', linestyles='--')
plt.clabel(C4, fontsize=6, colors='red')

# Contour lines for FDS 2 cm
C2 = plt.contour(X, Y, Z2, levels=levels, colors='green', linestyles='--')
plt.clabel(C2, fontsize=6, colors='green')

# Contour lines for FDS 1 cm
C1 = plt.contour(X, Y, Z1, levels=levels, colors='yellow', linestyles='--')
plt.clabel(C1, fontsize=6, colors='yellow')

# Add the legend
legend_elements = [mpatches.Patch(label='Experimental Data'), 
                   Line2D([0], [0], color='red', linestyle='--', label='FDS 4 cm'),
                   Line2D([0], [0], color='green', linestyle='--', label='FDS 2 cm'),
                   Line2D([0], [0], color='yellow', linestyle='--', label='FDS 1 cm')]
plt.legend(handles=legend_elements, loc='lower right', fontsize=8)

fig.savefig(pltdir + 'Marinite_Heat_Flux_Contours.pdf', format='pdf')
plt.close(fig)

