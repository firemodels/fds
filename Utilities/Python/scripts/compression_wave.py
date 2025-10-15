
# compression_wave case

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import fdsplotlib

plot_style = fdsplotlib.get_plot_style('fds')

outdir = '../../Verification/Scalar_Analytical_Solution/'
pltdir = '../../Manuals/FDS_Verification_Guide/SCRIPT_FIGURES/'

git_file = outdir + 'compression_wave_FL0_16_git.txt'
version_string = fdsplotlib.get_version_string(git_file)


def compression_wave_soln(rho0, x, y, a, c, t):

    b = np.sqrt(-1 + a**2)
    d = np.sqrt(-1 + c**2)

    x0 = 2 * np.arctan(b/a * np.tan(np.arctan((1 + a*np.tan(x/2))/b) - b*t/2) - 1/a)
    y0 = 2 * np.arctan(d/c * np.tan(np.arctan((1 + c*np.tan(y/2))/d) - d*t/2) - 1/c)

    Ix0 = np.log(-a**2 - np.cos(2*np.arctan((1 + a*np.tan(x0/2))/b)) + b*np.sin(2*np.arctan((1 + a*np.tan(x0/2))/b)))
    Iy0 = np.log(-c**2 - np.cos(2*np.arctan((1 + c*np.tan(y0/2))/d)) + d*np.sin(2*np.arctan((1 + c*np.tan(y0/2))/d)))

    Ix = np.log(-a**2 - np.cos(b*t + 2*np.arctan((1 + a*np.tan(x0/2))/b)) + b*np.sin(b*t + 2*np.arctan((1 + a*np.tan(x0/2))/b)))
    Iy = np.log(-c**2 - np.cos(d*t + 2*np.arctan((1 + c*np.tan(y0/2))/d)) + d*np.sin(d*t + 2*np.arctan((1 + c*np.tan(y0/2))/d)))

    q0 = np.log(rho0)
    q = q0 + Ix - Ix0 + Iy - Iy0

    rho = np.exp(q)

    return rho


# central differencing, FL=0
M_FL0_16 = pd.read_csv(os.path.join(outdir, 'compression_wave_FL0_16_devc.csv'), skiprows=2, header=None).values
M_FL0_32 = pd.read_csv(os.path.join(outdir, 'compression_wave_FL0_32_devc.csv'), skiprows=2, header=None).values
M_FL0_64 = pd.read_csv(os.path.join(outdir, 'compression_wave_FL0_64_devc.csv'), skiprows=2, header=None).values
M_FL0_128 = pd.read_csv(os.path.join(outdir, 'compression_wave_FL0_128_devc.csv'), skiprows=2, header=None).values
t_FL0_16 = M_FL0_16[:, 0]; rho_fds_FL0_16 = M_FL0_16[:, 2]
t_FL0_32 = M_FL0_32[:, 0]; rho_fds_FL0_32 = M_FL0_32[:, 2]
t_FL0_64 = M_FL0_64[:, 0]; rho_fds_FL0_64 = M_FL0_64[:, 2]
t_FL0_128 = M_FL0_128[:, 0]; rho_fds_FL0_128 = M_FL0_128[:, 2]

# Superbee limiter, FL=2
M_FL2_16 = pd.read_csv(os.path.join(outdir, 'compression_wave_FL2_16_devc.csv'), skiprows=2, header=None).values
M_FL2_32 = pd.read_csv(os.path.join(outdir, 'compression_wave_FL2_32_devc.csv'), skiprows=2, header=None).values
M_FL2_64 = pd.read_csv(os.path.join(outdir, 'compression_wave_FL2_64_devc.csv'), skiprows=2, header=None).values
M_FL2_128 = pd.read_csv(os.path.join(outdir, 'compression_wave_FL2_128_devc.csv'), skiprows=2, header=None).values
t_FL2_16 = M_FL2_16[:, 0]; rho_fds_FL2_16 = M_FL2_16[:, 2]
t_FL2_32 = M_FL2_32[:, 0]; rho_fds_FL2_32 = M_FL2_32[:, 2]
t_FL2_64 = M_FL2_64[:, 0]; rho_fds_FL2_64 = M_FL2_64[:, 2]
t_FL2_128 = M_FL2_128[:, 0]; rho_fds_FL2_128 = M_FL2_128[:, 2]

# CHARM limiter, FL=4
M_FL4_16 = pd.read_csv(os.path.join(outdir, 'compression_wave_FL4_16_devc.csv'), skiprows=2, header=None).values
M_FL4_32 = pd.read_csv(os.path.join(outdir, 'compression_wave_FL4_32_devc.csv'), skiprows=2, header=None).values
M_FL4_64 = pd.read_csv(os.path.join(outdir, 'compression_wave_FL4_64_devc.csv'), skiprows=2, header=None).values
M_FL4_128 = pd.read_csv(os.path.join(outdir, 'compression_wave_FL4_128_devc.csv'), skiprows=2, header=None).values
t_FL4_16 = M_FL4_16[:, 0]; rho_fds_FL4_16 = M_FL4_16[:, 2]
t_FL4_32 = M_FL4_32[:, 0]; rho_fds_FL4_32 = M_FL4_32[:, 2]
t_FL4_64 = M_FL4_64[:, 0]; rho_fds_FL4_64 = M_FL4_64[:, 2]
t_FL4_128 = M_FL4_128[:, 0]; rho_fds_FL4_128 = M_FL4_128[:, 2]

# analytical solution

a = 2
c = 3

L = 2*np.pi
x = 1.5*np.pi
y = 1.5*np.pi

rho_FL0_16 = compression_wave_soln(rho_fds_FL0_16[0], x-L/32, y-L/32, a, c, t_FL0_16)
error_FL0_16 = np.linalg.norm(rho_fds_FL0_16-rho_FL0_16)/len(t_FL0_16)
rho_FL0_32 = compression_wave_soln(rho_fds_FL0_32[0], x-L/64, y-L/64, a, c, t_FL0_32)
error_FL0_32 = np.linalg.norm(rho_fds_FL0_32-rho_FL0_32)/len(t_FL0_32)
rho_FL0_64 = compression_wave_soln(rho_fds_FL0_64[0], x-L/128, y-L/128, a, c, t_FL0_64)
error_FL0_64 = np.linalg.norm(rho_fds_FL0_64-rho_FL0_64)/len(t_FL0_64)
rho_FL0_128 = compression_wave_soln(rho_fds_FL0_128[0], x-L/256, y-L/256, a, c, t_FL0_128)
error_FL0_128 = np.linalg.norm(rho_fds_FL0_128-rho_FL0_128)/len(t_FL0_128)

rho_FL2_16 = compression_wave_soln(rho_fds_FL2_16[0], x-L/32, y-L/32, a, c, t_FL2_16)
error_FL2_16 = np.linalg.norm(rho_fds_FL2_16-rho_FL2_16)/len(t_FL2_16)
rho_FL2_32 = compression_wave_soln(rho_fds_FL2_32[0], x-L/64, y-L/64, a, c, t_FL2_32)
error_FL2_32 = np.linalg.norm(rho_fds_FL2_32-rho_FL2_32)/len(t_FL2_32)
rho_FL2_64 = compression_wave_soln(rho_fds_FL2_64[0], x-L/128, y-L/128, a, c, t_FL2_64)
error_FL2_64 = np.linalg.norm(rho_fds_FL2_64-rho_FL2_64)/len(t_FL2_64)
rho_FL2_128 = compression_wave_soln(rho_fds_FL2_128[0], x-L/256, y-L/256, a, c, t_FL2_128)
error_FL2_128 = np.linalg.norm(rho_fds_FL2_128-rho_FL2_128)/len(t_FL2_128)

rho_FL4_16 = compression_wave_soln(rho_fds_FL4_16[0], x-L/32, y-L/32, a, c, t_FL4_16)
error_FL4_16 = np.linalg.norm(rho_fds_FL4_16-rho_FL4_16)/len(t_FL4_16)
rho_FL4_32 = compression_wave_soln(rho_fds_FL4_32[0], x-L/64, y-L/64, a, c, t_FL4_32)
error_FL4_32 = np.linalg.norm(rho_fds_FL4_32-rho_FL4_32)/len(t_FL4_32)
rho_FL4_64 = compression_wave_soln(rho_fds_FL4_64[0], x-L/128, y-L/128, a, c, t_FL4_64)
error_FL4_64 = np.linalg.norm(rho_fds_FL4_64-rho_FL4_64)/len(t_FL4_64)
rho_FL4_128 = compression_wave_soln(rho_fds_FL4_128[0], x-L/256, y-L/256, a, c, t_FL4_128)
error_FL4_128 = np.linalg.norm(rho_fds_FL4_128-rho_FL4_128)/len(t_FL4_128)

fig = fdsplotlib.plot_to_fig(x_data=t_FL4_128, y_data=rho_FL4_128, marker_style='k-', data_label='Analytical Solution',
                             x_min=0, x_max=12.5, y_min=0, y_max=8,
                             revision_label=version_string,
                             x_label='Time (s)',
                             y_label='Density (kg/m$^3$)')

fdsplotlib.plot_to_fig(x_data=t_FL4_16, y_data=rho_fds_FL4_16, marker_style='c--', data_label='FDS $N=16$', figure_handle=fig)
fdsplotlib.plot_to_fig(x_data=t_FL4_32, y_data=rho_fds_FL4_32, marker_style='g--', data_label='FDS $N=32$', figure_handle=fig)
fdsplotlib.plot_to_fig(x_data=t_FL4_64, y_data=rho_fds_FL4_64, marker_style='b--', data_label='FDS $N=64$', figure_handle=fig)
fdsplotlib.plot_to_fig(x_data=t_FL4_128,y_data=rho_fds_FL4_128,marker_style='r--', data_label='FDS $N=128$',figure_handle=fig)

plt.savefig(os.path.join(pltdir, 'compression_wave_time_series.pdf'), format='pdf')
plt.close()

# convergence plot

h = 2*np.pi/np.array([16, 32, 64, 128])
e_FL0 = np.array([error_FL0_16, error_FL0_32, error_FL0_64, error_FL0_128])
e_FL2 = np.array([error_FL2_16, error_FL2_32, error_FL2_64, error_FL2_128])
e_FL4 = np.array([error_FL4_16, error_FL4_32, error_FL4_64, error_FL4_128])

fig = fdsplotlib.plot_to_fig(x_data=h, y_data=0.1*h, marker_style='k--', data_label=r'$O(\delta x)$',
                             x_min=1e-2, x_max=1, y_min=1e-4, y_max=1e-1,
                             plot_type='loglog',
                             revision_label=version_string,
                             x_label='Grid Spacing (m)',
                             y_label='L$_2$ Error (kg/m$^3$)')

fdsplotlib.plot_to_fig(x_data=h, y_data=0.1*h**2, marker_style='k-', data_label=r'$O(\delta x^2)$', figure_handle=fig)
fdsplotlib.plot_to_fig(x_data=h, y_data=e_FL0, marker_style='b*-', data_label='Central', figure_handle=fig)
fdsplotlib.plot_to_fig(x_data=h, y_data=e_FL2, marker_style='ro-', data_label='Superbee', figure_handle=fig)
fdsplotlib.plot_to_fig(x_data=h, y_data=e_FL4, marker_style='g^-', data_label='CHARM', figure_handle=fig)

plt.savefig(os.path.join(pltdir, 'compression_wave_convergence.pdf'), format='pdf')
plt.close()

