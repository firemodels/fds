
# Plot results of the Heskestad_Flame_Height cases against various experimental correlations

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
import matplotlib
import fdsplotlib

plot_style = fdsplotlib.get_plot_style('fds')

outdir = '../../../out/Heskestad_Flame_Height/'
expdir = '../../../exp/Heskestad_Flame_Height/'
pltdir = '../../Manuals/FDS_Validation_Guide/SCRIPT_FIGURES/Heskestad/'

git_file = outdir + 'Qs=p1_RI=05_git.txt'
version_string = fdsplotlib.get_version_string(git_file)

# list of line files
filename = [
    ['Qs=p1_RI=05_line.csv',    'Qs=p1_RI=10_line.csv',    'Qs=p1_RI=20_line.csv'],
    ['Qs=p2_RI=05_line.csv',    'Qs=p2_RI=10_line.csv',    'Qs=p2_RI=20_line.csv'],
    ['Qs=p5_RI=05_line.csv',    'Qs=p5_RI=10_line.csv',    'Qs=p5_RI=20_line.csv'],
    ['Qs=1_RI=05_line.csv',     'Qs=1_RI=10_line.csv',     'Qs=1_RI=20_line.csv'],
    ['Qs=2_RI=05_line.csv',     'Qs=2_RI=10_line.csv',     'Qs=2_RI=20_line.csv'],
    ['Qs=5_RI=05_line.csv',     'Qs=5_RI=10_line.csv',     'Qs=5_RI=20_line.csv'],
    ['Qs=10_RI=05_line.csv',    'Qs=10_RI=10_line.csv',    'Qs=10_RI=20_line.csv'],
    ['Qs=20_RI=05_line.csv',    'Qs=20_RI=10_line.csv',    'Qs=20_RI=20_line.csv'],
    ['Qs=50_RI=05_line.csv',    'Qs=50_RI=10_line.csv',    'Qs=50_RI=20_line.csv'],
    ['Qs=100_RI=05_line.csv',   'Qs=100_RI=10_line.csv',   'Qs=100_RI=20_line.csv'],
    ['Qs=200_RI=05_line.csv',   'Qs=200_RI=10_line.csv',   'Qs=200_RI=20_line.csv'],
    ['Qs=500_RI=05_line.csv',   'Qs=500_RI=10_line.csv',   'Qs=500_RI=20_line.csv'],
    ['Qs=1000_RI=05_line.csv',  'Qs=1000_RI=10_line.csv',  'Qs=1000_RI=20_line.csv'],
    ['Qs=2000_RI=05_line.csv',  'Qs=2000_RI=10_line.csv',  'Qs=2000_RI=20_line.csv'],
    ['Qs=5000_RI=05_line.csv',  'Qs=5000_RI=10_line.csv',  'Qs=5000_RI=20_line.csv'],
    ['Qs=10000_RI=05_line.csv', 'Qs=10000_RI=10_line.csv', 'Qs=10000_RI=20_line.csv']
]

rho_inf = 1.2
cp = 1
T_inf = 293
g = 9.81
D = 1.13
Qdot = [151, 303, 756, 1513, 3025, 7564, 15127, 30255, 75636, 151273, 302545, 756363, 1512725, 3025450, 7563625, 15127250]

# Initialize arrays
Qstar = np.zeros(16)
L_99 = np.zeros((16, 3))
L_95 = np.zeros((16, 3))

for i in range(16):  # hrr loop
    for j in range(3):  # resolution loop

        M = np.loadtxt(os.path.join(outdir, filename[i][j]), delimiter=',', skiprows=2)
        z = M[:, 0]
        dz = z[1] - z[0]
        hrrpul = M[:, 1]
        Qdot_line = np.sum(hrrpul) * dz
        Qstar[i] = Qdot[i] / (rho_inf * cp * T_inf * np.sqrt(g) * D**(5/2))

        # determine flame height
        hrr = np.zeros(len(z))
        for n in range(len(z)):
            hrr[n] = np.sum(hrrpul[:n+1]) * dz * Qdot[i] / Qdot_line  # cumulative heat release

        k_99 = np.where(hrr > (0.99) * Qdot[i])[0]
        if len(k_99) > 0:
            k = k_99[0]
            if k > 0:
                L_99[i, j] = z[k-1] + dz * ((0.99) * Qdot[i] - hrr[k-1]) / (hrr[k] - hrr[k-1])
            else:
                L_99[i, j] = dz * (0.99) * Qdot[i] / hrr[k]
        
        k_95 = np.where(hrr > (0.95) * Qdot[i])[0]
        if len(k_95) > 0:
            k = k_95[0]
            if k > 0:
                L_95[i, j] = z[k-1] + dz * ((0.95) * Qdot[i] - hrr[k-1]) / (hrr[k] - hrr[k-1])
            else:
                L_95[i, j] = dz * (0.95) * Qdot[i] / hrr[k]

# Set matplotlib configuration

fig = fdsplotlib.plot_to_fig(x_data=[1e-6,1e-6], y_data=[1e-6,1e-6],
                             x_min=0.05, x_max=2e4, y_min=0.1, y_max=3e2,
                             revision_label=version_string,
                             plot_title = 'Flame Height Variation',
                             plot_type='loglog',
                             x_label='$Q^*$',
                             y_label=r'$L_{\rm f}/D$')

# Load experimental data
M = pd.read_csv(os.path.join(expdir, 'flame_lengths.csv'), skiprows=0)
data = M.values

Steward = data[:, 4]
Becker_and_Liang = data[:, 5]
Cox_and_Chitty = data[:, 6]
Heskestad = data[:, 7]
Hasemi_and_Tokunaga = data[:, 8]
Cetegen = data[:, 9]
Delichatsios = data[:, 10]

fdsplotlib.plot_to_fig(x_data=Qstar, y_data=Steward, marker_style='k-', figure_handle=fig, data_label='Steward')
fdsplotlib.plot_to_fig(x_data=Qstar, y_data=Becker_and_Liang, marker_style='g-', figure_handle=fig, data_label=r'Becker & Liang')
fdsplotlib.plot_to_fig(x_data=Qstar[1:16], y_data=Cox_and_Chitty[1:16], marker_style='c-', figure_handle=fig, data_label=r'Cox & Chitty')
fdsplotlib.plot_to_fig(x_data=Qstar, y_data=Heskestad, marker_style='r-', figure_handle=fig, data_label='Heskestad')
fdsplotlib.plot_to_fig(x_data=Qstar[2:16], y_data=Hasemi_and_Tokunaga[2:16], marker_style='m-', figure_handle=fig, data_label=r'Hasemi & Tokunaga')
fdsplotlib.plot_to_fig(x_data=Qstar, y_data=Cetegen, marker_style='b-', figure_handle=fig, data_label='Cetegen')
fdsplotlib.plot_to_fig(x_data=Qstar, y_data=Delichatsios, marker_style='y-', figure_handle=fig, data_label='Delichatsios')
fdsplotlib.plot_to_fig(x_data=Qstar, y_data=np.max(L_99[:, :3], axis=1), marker_style='r--', figure_handle=fig, data_label='Max FDS 99 %')
fdsplotlib.plot_to_fig(x_data=Qstar, y_data=np.min(L_95[:, :3], axis=1), marker_style='b--', figure_handle=fig, data_label='Min FDS 99 %')

plt.savefig(pltdir + 'Flame_Height2.pdf', format='pdf')

