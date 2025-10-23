"""
McGrattan
12-14-2016
wind_profiles.py

Atmospheric Boundary Layer profiles, based on M-O theory as described in
"Falcon Series Data Report", GRI-89/0138, June 1990.
"""

import pandas as pd
import numpy as np
from fdsplotlib import plot_to_fig
from matplotlib import pyplot as plt
import os

# Directories
outdir = '../../../out/LNG_Dispersion/'
expdir = '../../../exp/LNG_Dispersion/'
pltdir = '../../Manuals/FDS_Validation_Guide/SCRIPT_FIGURES/LNG_Dispersion/'

# Physical constants
g = 9.81
cp = 1005
kappa = 0.41
p_0 = 1000

# Heights
z = np.array([0.001, 0.1, 0.5, 1, 1.9, 3, 4, 5, 6, 7, 8, 9, 10, 15, 20, 25, 30])

# Test-specific parameters
test = ['Burro3','Burro7','Burro8','Burro9',
        'Coyote3','Coyote5','Coyote6',
        'Falcon1','Falcon3','Falcon4',
        'MaplinSands27','MaplinSands34','MaplinSands35']

u_star =     [0.255, 0.372, 0.074, 0.252, 0.310, 0.480, 0.210, 0.0605, 0.3053, 0.3694, 0.19, 0.280, 0.315]
theta_star = [-0.532, -0.097, 0.026, -0.035, -0.890, -0.520, 0.039, 0.0577, -0.0175, 0.1521, -0.180, -0.0745, -0.0879]
L =          [-9.49, -110.8, 16.2, -142., -8.56, -33.2, 82.5, 4.963, -422.2, 69.38, -14.4, -75.5, -81.15]
p_r =        [948, 940, 941, 940., 936, 939, 942, 908.9, 900.8, 906.3, 1000., 1000., 1000.]
z_0 =        [0.0002]*7 + [0.008]*3 + [0.00002]*3
z_r =        [1,1,1,1,0.5,0.5,0.5,1,1,1,1.9,1.9,1.9]
i_z_r =      [4,4,4,4,3,3,3,4,4,4,5,5,5]

# Loop over tests
for i in range(len(test)):

    # Read CSVs and force numeric conversion
    M_path = f"{expdir}{test[i]}_profile.csv"
    S_path = f"{outdir}{test[i]}_line.csv"

    M = pd.read_csv(M_path, skiprows=2)
    S = pd.read_csv(S_path, skiprows=2)

    M = M.apply(pd.to_numeric, errors='coerce')
    S = S.apply(pd.to_numeric, errors='coerce')

    # Monin-Obukhov functions
    if L[i] < 0:
        x = (1 - 16*z/L[i])**0.25
        psi_m = 2*np.log((1+x)/2) + np.log((1+x**2)/2) - 2*np.arctan(x) + np.pi/2
        psi_h = 2*np.log((1+x**2)/2)
    else:
        psi_m = -5*z/L[i]
        psi_h = psi_m

    # Velocity and temperature calculations
    u = (u_star[i]/kappa)*(np.log(z/z_0[i]) - psi_m)
    T_a = M.iloc[:,2].values + 273.15
    theta_r = (p_0/p_r[i])**0.285*(T_a + (g/cp)*(M.iloc[:,0].values - z_r[i]))
    theta_0 = theta_r[0] - (theta_star[i]/kappa)*(np.log(z_r[i]/z_0[i]) - psi_h[i_z_r[i]-1])
    theta = theta_0 + (theta_star[i]/kappa)*(np.log(z/z_0[i]) - psi_h)
    T = theta*(p_0/p_r[i])**-0.285 - (g/cp)*(z - z_r[i])

    # Plot velocity
    fig_u = plot_to_fig(x_data=M.iloc[:,1].values, y_data=M.iloc[:,0].values, figure_handle=None, marker_style='ko',
                        x_label='Velocity (m/s)', y_label='Height (m)', plot_title=f'{test[i]} Velocity',
                        x_min=0,x_max=15,y_min=0,y_max=30, data_label='Measured', legend_location='lower right')
    plot_to_fig(x_data=u, y_data=z, figure_handle=fig_u, marker_style='k-', data_label='M-O Theory')
    plot_to_fig(x_data=S.iloc[:,1].values, y_data=S.iloc[:,0].values, figure_handle=fig_u, marker_style='k--', data_label='Simulated')

    fig_u.savefig(os.path.join(pltdir,f'{test[i]}_vel.pdf'))
    plt.close(fig_u)

    # Plot temperature
    col2_first_valid = M.iloc[:,2].dropna().values[0]  # first numeric value in column 2
    x_min = col2_first_valid - 5
    x_max = col2_first_valid + 5
    fig_T = plot_to_fig(x_data=M.iloc[:,2].values, y_data=M.iloc[:,0].values, figure_handle=None, marker_style='ko',
                        x_label='Temperature (°C)', y_label='Height (m)', plot_title=f'{test[i]} Temperature',
                        x_min=x_min, x_max=x_max, y_min=0, y_max=30, data_label='Measured', legend_location='lower right')
    plot_to_fig(x_data=T-273.15, y_data=z, figure_handle=fig_T, marker_style='k-', data_label='M-O Theory')
    plot_to_fig(x_data=S.iloc[:,2].values, y_data=S.iloc[:,0].values, figure_handle=fig_T, marker_style='k--', data_label='Simulated')

    fig_T.savefig(os.path.join(pltdir,f'{test[i]}_tmp.pdf'))
    plt.close(fig_T)

