
# Analytical Solution, Pulsating Wave

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import fdsplotlib

plot_style = fdsplotlib.get_plot_style('fds')

outdir = '../../Verification/Scalar_Analytical_Solution/'
pltdir = '../../Manuals/FDS_Verification_Guide/SCRIPT_FIGURES/'

git_file = outdir + 'pulsating_FL0_16_git.txt'
version_string = fdsplotlib.get_version_string(git_file)


def section2_soln(rho0, x, y, B, w, t):

    x0 = 2 * np.arctan(np.tan(x/2) * np.exp(-B/w * np.sin(w*t)))
    y0 = 2 * np.arctan(np.tan(y/2) * np.exp(-B/w * np.sin(w*t)))

    q0 = np.log(rho0)
    q = q0 + np.log((1 + (np.tan(x0/2))**2 * np.exp(2*B/w * np.sin(w*t))) / (1 + (np.tan(x0/2))**2)) \
           + np.log((1 + (np.tan(y0/2))**2 * np.exp(2*B/w * np.sin(w*t))) / (1 + (np.tan(y0/2))**2)) \
           - 2*B/w * np.sin(w*t)  # q(x,y,t) in Verification Guide

    rho = np.exp(q)

    return rho


devc_col = 3  # (x,y)=(pi,pi)-->devc_col=2, (x,y)=(1.5*pi,1.5*pi)-->devc_col=3

# Superbee limiter, FL=2
M_FL2_16 = np.genfromtxt(outdir + 'pulsating_FL2_16_devc.csv', delimiter=',', skip_header=2)
M_FL2_32 = np.genfromtxt(outdir + 'pulsating_FL2_32_devc.csv', delimiter=',', skip_header=2)
M_FL2_64 = np.genfromtxt(outdir + 'pulsating_FL2_64_devc.csv', delimiter=',', skip_header=2)
M_FL2_128 = np.genfromtxt(outdir + 'pulsating_FL2_128_devc.csv', delimiter=',', skip_header=2)
t_FL2_16 = M_FL2_16[:, 0]
rho_fds_FL2_16 = M_FL2_16[:, devc_col-1]
t_FL2_32 = M_FL2_32[:, 0]
rho_fds_FL2_32 = M_FL2_32[:, devc_col-1]
t_FL2_64 = M_FL2_64[:, 0]
rho_fds_FL2_64 = M_FL2_64[:, devc_col-1]
t_FL2_128 = M_FL2_128[:, 0]
rho_fds_FL2_128 = M_FL2_128[:, devc_col-1]

# CHARM limiter, FL=4
M_FL4_16 = np.genfromtxt(outdir + 'pulsating_FL4_16_devc.csv', delimiter=',', skip_header=2)
M_FL4_32 = np.genfromtxt(outdir + 'pulsating_FL4_32_devc.csv', delimiter=',', skip_header=2)
M_FL4_64 = np.genfromtxt(outdir + 'pulsating_FL4_64_devc.csv', delimiter=',', skip_header=2)
M_FL4_128 = np.genfromtxt(outdir + 'pulsating_FL4_128_devc.csv', delimiter=',', skip_header=2)
t_FL4_16 = M_FL4_16[:, 0]
rho_fds_FL4_16 = M_FL4_16[:, devc_col-1]
t_FL4_32 = M_FL4_32[:, 0]
rho_fds_FL4_32 = M_FL4_32[:, devc_col-1]
t_FL4_64 = M_FL4_64[:, 0]
rho_fds_FL4_64 = M_FL4_64[:, devc_col-1]
t_FL4_128 = M_FL4_128[:, 0]
rho_fds_FL4_128 = M_FL4_128[:, devc_col-1]

# central differencing, FL=0
M_FL0_16 = np.genfromtxt(outdir + 'pulsating_FL0_16_devc.csv', delimiter=',', skip_header=2)
M_FL0_32 = np.genfromtxt(outdir + 'pulsating_FL0_32_devc.csv', delimiter=',', skip_header=2)
M_FL0_64 = np.genfromtxt(outdir + 'pulsating_FL0_64_devc.csv', delimiter=',', skip_header=2)
M_FL0_128 = np.genfromtxt(outdir + 'pulsating_FL0_128_devc.csv', delimiter=',', skip_header=2)
t_FL0_16 = M_FL0_16[:, 0]
rho_fds_FL0_16 = M_FL0_16[:, devc_col-1]
t_FL0_32 = M_FL0_32[:, 0]
rho_fds_FL0_32 = M_FL0_32[:, devc_col-1]
t_FL0_64 = M_FL0_64[:, 0]
rho_fds_FL0_64 = M_FL0_64[:, devc_col-1]
t_FL0_128 = M_FL0_128[:, 0]
rho_fds_FL0_128 = M_FL0_128[:, devc_col-1]

# analytical solution

B = 1
w = 1

L = 2*np.pi
if devc_col == 2:
    x = np.pi
    y = np.pi
if devc_col == 3:
    x = 1.5*np.pi
    y = 1.5*np.pi

rho_FL2_16 = section2_soln(rho_fds_FL2_16[0], x-L/32, y-L/32, B, w, t_FL2_16)
error_FL2_16 = np.linalg.norm(rho_fds_FL2_16 - rho_FL2_16) / len(t_FL2_16)
rho_FL2_32 = section2_soln(rho_fds_FL2_32[0], x-L/64, y-L/64, B, w, t_FL2_32)
error_FL2_32 = np.linalg.norm(rho_fds_FL2_32 - rho_FL2_32) / len(t_FL2_32)
rho_FL2_64 = section2_soln(rho_fds_FL2_64[0], x-L/128, y-L/128, B, w, t_FL2_64)
error_FL2_64 = np.linalg.norm(rho_fds_FL2_64 - rho_FL2_64) / len(t_FL2_64)
rho_FL2_128 = section2_soln(rho_fds_FL2_128[0], x-L/256, y-L/256, B, w, t_FL2_128)
error_FL2_128 = np.linalg.norm(rho_fds_FL2_128 - rho_FL2_128) / len(t_FL2_128)

rho_FL4_16 = section2_soln(rho_fds_FL4_16[0], x-L/32, y-L/32, B, w, t_FL4_16)
error_FL4_16 = np.linalg.norm(rho_fds_FL4_16 - rho_FL4_16) / len(t_FL4_16)
rho_FL4_32 = section2_soln(rho_fds_FL4_32[0], x-L/64, y-L/64, B, w, t_FL4_32)
error_FL4_32 = np.linalg.norm(rho_fds_FL4_32 - rho_FL4_32) / len(t_FL4_32)
rho_FL4_64 = section2_soln(rho_fds_FL4_64[0], x-L/128, y-L/128, B, w, t_FL4_64)
error_FL4_64 = np.linalg.norm(rho_fds_FL4_64 - rho_FL4_64) / len(t_FL4_64)
rho_FL4_128 = section2_soln(rho_fds_FL4_128[0], x-L/256, y-L/256, B, w, t_FL4_128)
error_FL4_128 = np.linalg.norm(rho_fds_FL4_128 - rho_FL4_128) / len(t_FL4_128)

rho_FL0_16 = section2_soln(rho_fds_FL0_16[0], x-L/32, y-L/32, B, w, t_FL0_16)
error_FL0_16 = np.linalg.norm(rho_fds_FL0_16 - rho_FL0_16) / len(t_FL0_16)
rho_FL0_32 = section2_soln(rho_fds_FL0_32[0], x-L/64, y-L/64, B, w, t_FL0_32)
error_FL0_32 = np.linalg.norm(rho_fds_FL0_32 - rho_FL0_32) / len(t_FL0_32)
rho_FL0_64 = section2_soln(rho_fds_FL0_64[0], x-L/128, y-L/128, B, w, t_FL0_64)
error_FL0_64 = np.linalg.norm(rho_fds_FL0_64 - rho_FL0_64) / len(t_FL0_64)
rho_FL0_128 = section2_soln(rho_fds_FL0_128[0], x-L/256, y-L/256, B, w, t_FL0_128)
error_FL0_128 = np.linalg.norm(rho_fds_FL0_128 - rho_FL0_128) / len(t_FL0_128)

fig = fdsplotlib.plot_to_fig(x_data=t_FL2_128, y_data=rho_FL2_128, marker_style='k-', data_label='Analytical Solution',
                             x_min=0, x_max=12.5, y_min=0, y_max=2,
                             revision_label=version_string,
                             x_label='Time (s)',
                             y_label='Density (kg/m$^3$)')

fdsplotlib.plot_to_fig(x_data=t_FL2_16, y_data=rho_fds_FL2_16, marker_style='c--', data_label='FDS $N=16$', figure_handle=fig)
fdsplotlib.plot_to_fig(x_data=t_FL2_32, y_data=rho_fds_FL2_32, marker_style='g--', data_label='FDS $N=32$', figure_handle=fig)
fdsplotlib.plot_to_fig(x_data=t_FL2_64, y_data=rho_fds_FL2_64, marker_style='b--', data_label='FDS $N=64$', figure_handle=fig)
fdsplotlib.plot_to_fig(x_data=t_FL2_128,y_data=rho_fds_FL2_128,marker_style='r--', data_label='FDS $N=128$',figure_handle=fig)

plt.savefig(pltdir + 'pulsating_time_series.pdf', format='pdf')
plt.close()

# convergence plot

h = 2*np.pi / np.array([16, 32, 64, 128])

fig = fdsplotlib.plot_to_fig(x_data=h, y_data=0.1*h, marker_style='k--', data_label=r'$O(\delta x)$',
                             x_min=1e-2, x_max=1, y_min=1e-6, y_max=1e-1,
                             plot_type='loglog',
                             revision_label=version_string,
                             x_label='Grid Spacing (m)',
                             y_label='L$_2$ Error (kg/m$^3$)')

e_FL2 = np.array([error_FL2_16, error_FL2_32, error_FL2_64, error_FL2_128])
e_FL4 = np.array([error_FL4_16, error_FL4_32, error_FL4_64, error_FL4_128])
e_FL0 = np.array([error_FL0_16, error_FL0_32, error_FL0_64, error_FL0_128])

fdsplotlib.plot_to_fig(x_data=h, y_data=0.1*h**2, marker_style='k-', data_label=r'$O(\delta x^2)$', figure_handle=fig)
fdsplotlib.plot_to_fig(x_data=h, y_data=e_FL0, marker_style='b*-', data_label='FDS Central', figure_handle=fig)
fdsplotlib.plot_to_fig(x_data=h, y_data=e_FL2, marker_style='ro-', data_label='FDS Superbee', figure_handle=fig)
fdsplotlib.plot_to_fig(x_data=h, y_data=e_FL4, marker_style='g^-', data_label='FDS CHARM', figure_handle=fig)

plt.savefig(pltdir + 'pulsating_convergence.pdf', format='pdf')
plt.close()

