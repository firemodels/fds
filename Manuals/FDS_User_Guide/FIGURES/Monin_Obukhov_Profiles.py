
# Atmospheric Boundary Layer profiles, based on M-O theory as described in
# "Falcon Series Data Report", GRI-89/0138, June 1990.
# These plots are used as illustrations in the FDS User's Guide.

import numpy as np
import matplotlib.pyplot as plt
import fdsplotlib

plot_style = fdsplotlib.get_plot_style('fds')

rho_0 = 1.2
g = 9.81
cp = 1005
kappa = 0.41
p_0 = 100000

L   = np.array([-100,100,-350])
Ls  = np.array(['-100','100','-350'])
z_0 = np.array([0.01,0.005,0.005])
z_r = 2.0
u_r = 5.0

for i in range(3):

    p_r = p_0 - rho_0 * g * (z_r - z_0[i])
    T_r = 20 + 273.15
    theta_r = T_r * (p_0 / p_r) ** 0.285
    
    if L[i] < 0:
        x_r = (1 - 16 * z_r / L[i]) ** 0.25
        psi_m_r = 2 * np.log((1 + x_r) / 2) + np.log((1 + x_r**2) / 2) - 2 * np.arctan(x_r) + np.pi / 2
        psi_h_r = 2 * np.log((1 + x_r**2) / 2)
    else:
        psi_m_r = -5 * z_r / L[i]
        psi_h_r = psi_m_r
    
    u_star = kappa * u_r / (np.log(z_r / z_0[i]) - psi_m_r)
    theta_0 = theta_r / (1 + u_star**2 * (np.log(z_r / z_0[i]) - psi_h_r) / (g * kappa**2 * L[i]))
    theta_star = u_star**2 * theta_0 / (g * kappa * L[i])
    
    z = np.array([z_0[i], 10*z_0[i], 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 15, 20, 25, 30, 50, 100])
    
    # Figure 1 - Velocity profile
    
    if L[i] < 0:
        x = (1 - 16 * z / L[i]) ** 0.25
        psi_m = 2 * np.log((1 + x) / 2) + np.log((1 + x**2) / 2) - 2 * np.arctan(x) + np.pi / 2
        psi_h = 2 * np.log((1 + x**2) / 2)
    else:
        psi_m = -5 * z / L[i]
        psi_h = psi_m
    
    u = (u_star / kappa) * (np.log(z / z_0[i]) - psi_m)
    theta = theta_0 + (theta_star / kappa) * (np.log(z / z_0[i]) - psi_h)
    T = theta * (p_0 / (p_0 - rho_0 * g * (z - z_0[i]))) ** -0.285
    
    fig = fdsplotlib.plot_to_fig(x_data=u, y_data=z, marker_style='ko-', 
                                 x_label='Velocity (m/s)', y_label='Height (m)',
                                 x_min=0, x_max=10, y_min=0, y_max=30)
    
    ax1 = fig.axes[0]
    ax1.text(0.05, 0.90, r'$u(2 \; \mathrm{m})=5$ m/s', transform=ax1.transAxes, verticalalignment='top')
    ax1.text(0.05, 0.80, f'$L={L[i]:4.0f}$ m', transform=ax1.transAxes, verticalalignment='top')
    ax1.text(0.05, 0.70, f'$z_0={z_0[i]:5.3f}$ m', transform=ax1.transAxes, verticalalignment='top')
    
    plt.savefig(f'Monin_Obukhov_vel_L_{Ls[i]}.pdf', format='pdf')
    plt.close()
    
    # Figure 2 - Temperature profile
    
    fig = fdsplotlib.plot_to_fig(x_data=T-273.15, y_data=z, marker_style='ko-',
                                 x_label='Temperature (°C)', y_label='Height (m)',
                                 x_min=15, x_max=25, y_min=0, y_max=30)
    
    ax2 = fig.axes[0]
    ax2.text(0.05, 0.90, r'$T(2 \; \mathrm{m})=20\;^\circ$C', transform=ax2.transAxes, verticalalignment='top')
    ax2.text(0.05, 0.80, f'$L={L[i]:4.0f}$ m', transform=ax2.transAxes, verticalalignment='top')
    ax2.text(0.05, 0.70, f'$z_0={z_0[i]:5.3f}$ m', transform=ax2.transAxes, verticalalignment='top')
    
    plt.savefig(f'Monin_Obukhov_tmp_L_{Ls[i]}.pdf', format='pdf')
    plt.close()

