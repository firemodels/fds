
# Atmospheric Boundary Layer profiles, based on M-O theory as described in
# "Falcon Series Data Report", GRI-89/0138, June 1990.
# These plots are used as illustrations in the FDS User's Guide.

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import fdsplotlib

plot_style = fdsplotlib.get_plot_style('fds')

outdir = '../../Verification/Atmospheric_Effects/'
pltdir = '../../Manuals/FDS_User_Guide/SCRIPT_FIGURES/'

git_file = outdir + 'atmospheric_boundary_layer_1_git.txt'
version_string = fdsplotlib.get_version_string(git_file)

basein  = outdir + 'atmospheric_boundary_layer'
baseout = pltdir + 'atmospheric_boundary_layer'

for i in range(1, 5):  # Loop from 1 to 4 (inclusive)

    datafile = basein + '_' + str(i)
    outfile = baseout + '_' + str(i)

    M1 = pd.read_csv(datafile + '_devc.csv', sep=',', skiprows=1)
    M2 = pd.read_csv(datafile + '_line.csv', sep=',', skiprows=1)

    u_r = M1.iloc[-1, 1]
    T_r = M1.iloc[-1, 2] + 273.15

    rho_0 = 1.2
    g = 9.81
    cp = 1005
    kappa = 0.41
    p_0 = 100000
    qdot = {1: 50, 2: -50, 3: 25, 4: -25}
    z_0 = {1: 0.25, 2: 0.25, 3: 0.125, 4: 0.125}
    T_low = {1: 15, 2: 15, 3: 15, 4: 15}
    T_high = {1: 25, 2: 25, 3: 25, 4: 25}
    u_high = {1: 20, 2: 20, 3: 10, 4: 15}
    fvec = {1: 0.01, 2: 0.01, 3: 0.002, 4: 0.005}
    s = {1: 8.15, 2: 8.15, 3: 4.075, 4: 4.075}

    theta_0 = T_r
    z_r = 10.
    p_r = p_0 - rho_0 * g * (z_r - z_0[i])
    theta_r = T_r * (p_0 / p_r) ** 0.285
    u_star = kappa * u_r / np.log(z_r / z_0[i])
    L = -u_star**3 * theta_0 * rho_0 * cp / (g * kappa * qdot[i])
    theta_star = u_star**2 * theta_0 / (g * kappa * L)

    z = np.array([z_0[i], 10*z_0[i], 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 15, 20, 25, 30, 50, 100])

    # Create figure 1 for velocity

    if L < 0:
        x = (1 - 16*z/L)**0.25
        psi_m = 2*np.log((1+x)/2) + np.log((1+x**2)/2) - 2*np.arctan(x) + np.pi/2
        psi_h = 2*np.log((1+x**2)/2)
    else:
        psi_m = -5*z/L
        psi_h = psi_m

    u = (u_star/kappa) * (np.log(z/z_0[i]) - psi_m)
    theta = theta_0 + (theta_star/kappa) * (np.log(z/z_0[i]) - psi_h)
    T = theta * (p_0 / (p_0 - rho_0*g*(z - z_0[i])))**(-0.285)

    T = T + (theta_0 - T[11]) 

    ERROR = abs(u[-1] - M2.iloc[-1, 1])
    if ERROR > 2.:
        print(f'Python Warning: atmospheric_boundary_layer Case {i} velocity out of tolerance. ERROR = {ERROR} m/s')

    fig = fdsplotlib.plot_to_fig(x_data=M2.iloc[:, 1].values, y_data=M2.iloc[:, 0].values, marker_style='k-', data_label='FDS',
                                 x_label='Velocity (m/s)', y_label='Height (m)',
                                 x_min=0, x_max=u_high[i], y_min=0, y_max=100)

    fdsplotlib.plot_to_fig(x_data=u, y_data=z, figure_handle=fig, marker_style='ko', data_label='M-O Theory')

    ax1 = fig.axes[0]

    ax1.text(0.05, 0.90, f'Case {i}', transform=ax1.transAxes)
    ax1.text(0.05, 0.80, f'$F={fvec[i]:.3f}' + r'\; \mathrm{Pa/m}$', transform=ax1.transAxes)
    ax1.text(0.05, 0.70, f'$s={s[i]:.2f}' + r'\; \mathrm{m}$', transform=ax1.transAxes)
    ax1.text(0.05, 0.60, f'$\\dot{{q}}_{{\\rm c}}\\prime \\prime={qdot[i]/1000:.3f}' + r'\; \mathrm{kW/m}^2$', transform=ax1.transAxes)
    ax1.text(0.05, 0.50, f'$u({z_r:.0f}' + r'\; \mathrm{m})=' + f'{u_r:.1f}' + r'\; \mathrm{m/s}$', transform=ax1.transAxes)
    ax1.text(0.05, 0.40, f'$L={L:.0f}$ m', transform=ax1.transAxes)
    ax1.text(0.05, 0.30, f'$z_0={z_0[i]:.3f}$ m', transform=ax1.transAxes)

    plt.savefig(outfile + '_vel.pdf', format='pdf')
    plt.close()

    # Create figure 2 for temperature

    ERROR = abs(T[-1] - 273.15 - M2.iloc[-1, 2])
    if ERROR > 1.0:
        print(f'Python Warning: atmospheric_boundary_layer Case {i} temperature out of tolerance. ERROR = {ERROR} K')

    fig = fdsplotlib.plot_to_fig(x_data=M2.iloc[:,2].values, y_data=M2.iloc[:,0].values, marker_style='k-', data_label='FDS',
                                 x_label='Temperature (°C)', y_label='Height (m)', 
                                 x_min=T_low[i], x_max=T_high[i], y_min=0, y_max=100)

    fdsplotlib.plot_to_fig(x_data=T - 273.15, y_data=z, figure_handle=fig, marker_style='ko', data_label='M-O Theory')

    ax2 = fig.axes[0]

    ax2.text(0.05, 0.90, f'Case {i}', transform=ax2.transAxes)
    ax2.text(0.05, 0.80, f'$T({z_r:.0f}' + r'\; \mathrm{m})=' + f'{T_r-273:.1f}' + r'\;^\circ$C', transform=ax2.transAxes)

    plt.savefig(outfile + '_tmp.pdf', format='pdf')
    plt.close()

