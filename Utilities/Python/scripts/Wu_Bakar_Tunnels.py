
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
import warnings
import fdsplotlib
warnings.filterwarnings('ignore')

# McGrattan
# 07-02-2025  
# Wu_Bakar_Tunnels.py
#
# This script replicates Fig. 4 of Wu and Bakar, Fire Safety Journal, 35 (2000) 363-390.

# Clear all variables and close all figures (Python equivalent)
plt.close('all')

outdir = '../../../out/Wu_Bakar_Tunnels/'
pltdir = '../../Manuals/FDS_Validation_Guide/SCRIPT_FIGURES/Wu_Bakar_Tunnels/'

g = 9.81
cp = 1.0
Tinf = 293
rho = 1.2

label = ['A', 'B', 'C', 'D', 'E']  # Tunnel names

Hbar = np.array([0.176, 0.250, 0.333, 0.400, 0.238])  # Hydraulic tunnel diameters, 4*A/P
Q = np.array([1.5, 3.0, 7.5, 10.5, 12.0, 15.0, 22.5, 30.0])  # HRR (kW) of the propane burner

V = np.array([[0.43, 0.46, 0.48, 0.48, 0.48, 0.48, 0.48, 0.48],
              [0.39, 0.48, 0.56, 0.59, 0.60, 0.60, 0.60, 0.60],
              [0.37, 0.45, 0.54, 0.57, 0.59, 0.60, 0.62, 0.65],
              [0.34, 0.40, 0.50, 0.54, 0.56, 0.59, 0.65, 0.65],
              [0.44, 0.54, 0.60, 0.60, 0.60, 0.60, 0.60, 0.60]])  # Measured critical velocities (m/s)

L = np.zeros((5, 8))

for i in range(5):
    # Find the back-layer length, L, defined as the distance from the upstream edge of the burner, 6.15 m, back to where
    # the back-layer temperature drops to 30 C, or 10 C above ambient.
    
    try:
        S = pd.read_csv(outdir + 'Wu_Bakar_Tunnel_' + label[i] + '_line.csv', skiprows=2)
        data = S.values
        
        for j in range(8):
            indices = np.where((data[:, j+1] > 30) & (data[:, j+1] < 50))[0]
            if len(indices) > 0:
                L[i, j] = 6.15 - data[indices[0], 0]
    except:
        print(f"Warning: Could not load data for Tunnel {label[i]}")

# Estimate the simulated critical velocity, Vmod, based on Li, Lei, and Ingason, FSJ 45 (2010) 361-370, Eq. 30.

Vmod = np.zeros((5, 8))
Qstar = np.zeros((5, 8))
Vstar = np.zeros((5, 8))

for i in range(5):
    for j in range(8):
        Vmod[i, j] = V[i, j] * (1 + L[i, j] / (18.5 * Hbar[i]))
        Qstar[i, j] = Q[j] / (rho * cp * Tinf * np.sqrt(g) * Hbar[i]**2.5)
        Vstar[i, j] = Vmod[i, j] / np.sqrt(g * Hbar[i])

# Difference between experimental critical velocity, V, and predicted, Vmod

Difference = 100 * (Vmod - V) / V

# Correlation of Wu and Bakar, Eq. 15-16, FSJ, 2000

Qcorr = np.arange(0.01, 2.05, 0.05)
Vcorr = np.zeros(len(Qcorr))

for i in range(len(Qcorr)):
    if Qcorr[i] < 0.2:
        Vcorr[i] = 0.4 * 0.2**(-1/3) * Qcorr[i]**(1/3)
    if Qcorr[i] >= 0.2:
        Vcorr[i] = 0.4

# Plot Qstar vs Vstar along with Qcorr vs Vcorr (Wu and Bakar, Fig. 4)

plot_style = fdsplotlib.get_plot_style('fds')
git_file = outdir + 'Wu_Bakar_Tunnel_A_git.txt'
version_string = fdsplotlib.get_version_string(git_file)

fig = fdsplotlib.plot_to_fig(x_data=Qcorr, y_data=Vcorr, data_label='Correlation', marker_style='k-',
                             x_min=0.001, x_max=10, y_min=0.1, y_max=1,
                             revision_label=version_string,
                             x_label='Heat Release Rate, $Q^*$', y_label='Critical Velocity, $V^*$',
                             legend_location='upper left',
                             plot_type='loglog')

fdsplotlib.plot_to_fig(x_data=Qstar[0, :], y_data=Vstar[0, :], marker_style='kd', data_label='Tunnel A', figure_handle=fig)
fdsplotlib.plot_to_fig(x_data=Qstar[1, :], y_data=Vstar[1, :], marker_style='rs', data_label='Tunnel B', figure_handle=fig)
fdsplotlib.plot_to_fig(x_data=Qstar[2, :], y_data=Vstar[2, :], marker_style='m^', data_label='Tunnel C', figure_handle=fig)
fdsplotlib.plot_to_fig(x_data=Qstar[3, :], y_data=Vstar[3, :], marker_style='c+', data_label='Tunnel D', figure_handle=fig)
fdsplotlib.plot_to_fig(x_data=Qstar[4, :], y_data=Vstar[4, :], marker_style='go', data_label='Tunnel E', figure_handle=fig)

# Make the pdf figure

plt.savefig(pltdir + 'Wu_Bakar_Critical_Velocity.pdf', format='pdf')

# Clear Qstar and Vstar for experimental data

Qstar = np.zeros((5, 8))
Vstar = np.zeros((5, 8))

Qstar[0, :] = np.array([0.108, 0.216, 0.539, 0.754, 0.862, 1.078, 1.616, 2.155])
Qstar[1, :] = np.array([0.041, 0.082, 0.205, 0.288, 0.329, 0.411, 0.616, 0.822])
Qstar[2, :] = np.array([0.021, 0.041, 0.103, 0.144, 0.164, 0.205, 0.308, 0.41])
Qstar[3, :] = np.array([0.01, 0.02, 0.05, 0.07, 0.08, 0.1, 0.15, 1000])
Qstar[4, :] = np.array([0.046, 0.093, 0.232, 0.325, 0.372, 0.465, 0.67, 0.929])

Vstar[0, :] = np.array([0.333, 0.356, 0.372, 0.372, 0.372, 0.372, 0.372, 0.372])
Vstar[1, :] = np.array([0.248, 0.309, 0.354, 0.379, 0.382, 0.382, 0.382, 0.382])
Vstar[2, :] = np.array([0.203, 0.251, 0.298, 0.319, 0.325, 0.332, 0.344, 0.359])
Vstar[3, :] = np.array([0.162, 0.192, 0.241, 0.259, 0.27, 0.284, 0.313, 1000])
Vstar[4, :] = np.array([0.288, 0.353, 0.393, 0.393, 0.393, 0.393, 0.393, 0.393])

# Plot Wu and Bakar data along with Qcorr vs Vcorr (Wu and Bakar, Fig. 4)

fig2 = fdsplotlib.plot_to_fig(x_data=Qcorr, y_data=Vcorr, data_label='Correlation', marker_style='k-',
                              x_min=0.001, x_max=10, y_min=0.1, y_max=1,
                              revision_label=version_string,
                              x_label='Heat Release Rate, $Q^*$', y_label='Critical Velocity, $V^*$',
                              legend_location='upper left',
                              plot_type='loglog')

fdsplotlib.plot_to_fig(x_data=Qstar[0, :], y_data=Vstar[0, :], marker_style='kd', data_label='Tunnel A', figure_handle=fig2)
fdsplotlib.plot_to_fig(x_data=Qstar[1, :], y_data=Vstar[1, :], marker_style='rs', data_label='Tunnel B', figure_handle=fig2)
fdsplotlib.plot_to_fig(x_data=Qstar[2, :], y_data=Vstar[2, :], marker_style='m^', data_label='Tunnel C', figure_handle=fig2)
fdsplotlib.plot_to_fig(x_data=Qstar[3, :], y_data=Vstar[3, :], marker_style='c+', data_label='Tunnel D', figure_handle=fig2)
fdsplotlib.plot_to_fig(x_data=Qstar[4, :], y_data=Vstar[4, :], marker_style='go', data_label='Tunnel E', figure_handle=fig2)

# Make the pdf figure
plt.savefig(pltdir + 'Wu_Bakar_Critical_Velocity_Exp_Data.pdf', format='pdf')



