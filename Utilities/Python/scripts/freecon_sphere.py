
# FDS Verification Guide, free convection from a sphere (free_conv_sphere.pdf)

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import fdsplotlib

plot_style = fdsplotlib.get_plot_style('fds')

outdir = '../../../out/Convection/'
pltdir = '../../Manuals/FDS_Verification_Guide/SCRIPT_FIGURES/'

git_file = outdir + 'free_conv_sphere_1_16_git.txt'
version_string = fdsplotlib.get_version_string(git_file)

# calculations below were used for input file setup

g = 9.80665
T1 = 303
T2 = 293
Tm = 0.5 * (T1 + T2)
beta = 1.0 / Tm
MW = 28.85476  # FDS 'LJ AIR'
P0 = 101325
rho = P0 * MW / (8341.5 * Tm)
mu = 1.8216e-5
cp = 1000
k = 0.018216  # for Pr=1 fluid

Pr = cp * mu / k
nu = mu / rho
alpha = k / (rho * cp)

setup = 1

if setup:

    # see J.P. Holman p. 357 for correlations

    d1 = 0.001
    Ra1 = (g * beta * (T1 - T2) * d1**3) / (alpha * nu)
    Nu1 = 2 + 0.43 * Ra1**(0.25)
    Tau1 = 1.0 / Nu1 * d1**2 / alpha

    d2 = 0.01
    Ra2 = (g * beta * (T1 - T2) * d2**3) / (alpha * nu)
    Nu2 = 2 + 0.43 * Ra1**(0.25)
    Tau2 = 1.0 / Nu2 * d2**2 / alpha

    d3 = 0.1
    Ra3 = (g * beta * (T1 - T2) * d3**3) / (alpha * nu)
    Nu3 = 2 + 0.43 * Ra3**(0.25)
    Tau3 = 1.0 / Nu3 * d3**2 / alpha

    d4 = 1
    Ra4 = (g * beta * (T1 - T2) * d4**3) / (alpha * nu)
    Nu4 = 2 + 0.5 * Ra4**(0.25)
    Tau4 = 1.0 / Nu4 * d4**2 / alpha

    RAYLEIGH_1 = np.logspace(0, 5, 100)   # Yuge
    RAYLEIGH_2 = np.logspace(5, 10, 100)  # Amato and Tien

    NUSSELT_1 = 2.0 + 0.43 * RAYLEIGH_1**(0.25)
    NUSSELT_2 = 2.0 + 0.50 * RAYLEIGH_2**(0.25)

    fig = fdsplotlib.plot_to_fig(x_data=RAYLEIGH_1, y_data=NUSSELT_1, data_label='Yuge (1960)', line_style='k-',
                                 x_label='Rayleigh Number', y_label='Nusselt Number',
                                 x_min=1, x_max=1e10, y_min=1, y_max=5e2,
                                 revision_label=version_string,
                                 plot_type='loglog')

    fdsplotlib.plot_to_fig(x_data=RAYLEIGH_2, y_data=NUSSELT_2, figure_handle=fig, data_label='Amato & Tien (1972)', line_style='k--')

casename = [
    'free_conv_sphere_1',
    'free_conv_sphere_2',
    'free_conv_sphere_3',
    'free_conv_sphere_4'
]

delta = np.array([0.001, 0.01, 0.1, 1])
T = np.array([303, 303, 303, 303])
A = 4 * np.pi * (delta / 2)**2

marker_style = ['r^', 'bsq']
res = ['8', '16']

Ra = np.zeros(len(delta))

for j in range(len(res)):
    for i in range(len(delta)):

        M = pd.read_csv(outdir + casename[i] + '_' + res[j] + '_devc.csv', skiprows=1)
        t = M['Time'].values
        Q = np.mean(M['Q'].values[round(len(M['Q'].values) / 2):]) * 1000
        alpha = k / (rho * cp)
        nu = mu / rho
        b = 2.0 / (T[i] + T2)
        Ra[i] = (g * b * (T[i] - T2) * delta[i]**3) / (alpha * nu)
        Nu_FDS = (Q / A[i]) * (delta[i] / k) / (T2 - T[i])

        fdsplotlib.plot_to_fig(x_data=np.array([Ra[i]]), y_data=np.array([Nu_FDS]), figure_handle=fig,
                        data_label=f'FDS $\\mathit{{D/\\delta x}}$={res[j]}' if j == 0 and i == 0 else (f'FDS $\\mathit{{D/\\delta x}}$={res[j]}' if i == 0 else None),
                        marker_style=marker_style[j])

plt.savefig(pltdir + 'free_conv_sphere.pdf', format='pdf')
plt.close()

