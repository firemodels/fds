
# Synthetic Eddy Method

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import fdsplotlib

plot_style = fdsplotlib.get_plot_style('fds')

outdir = '../../Verification/Turbulence/'
pltdir = '../../Manuals/FDS_User_Guide/SCRIPT_FIGURES/'

git_file = outdir + 'sem_flat_leddy_p2_git.txt'
version_string = fdsplotlib.get_version_string(git_file)


error_tolerance = 0.01

# Flat profile

M = pd.read_csv(outdir + 'sem_flat_leddy_p2_line.csv', header=1)

k = M.columns.get_loc('umean')
z = M.iloc[:, k-1].values
umean = M.iloc[:, k].values

u0 = 1
uprof = u0 * np.ones(len(z))

fig = fdsplotlib.plot_to_fig(x_data=uprof, y_data=z, marker_style='k-', data_label='Prescribed mean',
                             revision_label=version_string,
                             x_label='$u$ (m/s)', y_label='$z$ (m)',
                             x_min=0.4, x_max=1.2, y_min=0, y_max=1)

fdsplotlib.plot_to_fig(x_data=1.1*uprof, y_data=z, figure_handle=fig, marker_style='k--', data_label='Prescribed rms')
fdsplotlib.plot_to_fig(x_data=0.9*uprof, y_data=z, figure_handle=fig, marker_style='k--', data_label=None)
fdsplotlib.plot_to_fig(x_data=umean, y_data=z, figure_handle=fig, marker_style='b>-', data_label='FDS mean')

k = M.columns.get_loc('urms')
z = M.iloc[:, k-1].values
urms = M.iloc[:, k].values

fdsplotlib.plot_to_fig(x_data=umean+urms, y_data=z, figure_handle=fig, marker_style='r--', data_label='FDS rms')
fdsplotlib.plot_to_fig(x_data=umean-urms, y_data=z, figure_handle=fig, marker_style='r--', data_label=None)

plt.savefig(pltdir + 'sem_flat_leddy_p2.pdf', format='pdf')
plt.close()

umean_error = np.linalg.norm(umean - uprof) / u0 / len(umean)
if umean_error > error_tolerance:
    print(f'Matlab Warning: sem_flat_leddy_p2.fds umean_error = {umean_error}')

urms_error = np.linalg.norm(umean + urms - 1.1*uprof) / u0 / len(umean)
if urms_error > error_tolerance:
    print(f'Matlab Warning: sem_flat_leddy_p2.fds urms_error = {urms_error}')


# Parabolic profile

M = pd.read_csv(outdir + 'sem_par_leddy_p2_line.csv', header=1)

k = M.columns.get_loc('umean')
z = M.iloc[:, k-1].values
umean = M.iloc[:, k].values

umax = 1
h = 0.5
a = umax / h**2
uprof = umax - a * (z - h)**2

fig = fdsplotlib.plot_to_fig(x_data=uprof, y_data=z, marker_style='k-', data_label='Prescribed mean',
                             revision_label=version_string,
                             x_label='$u$ (m/s)', y_label='$z$ (m)',
                             x_min=0, x_max=1.2, y_min=0, y_max=1)

fdsplotlib.plot_to_fig(x_data=1.1*uprof, y_data=z, figure_handle=fig, marker_style='k--', data_label='Prescribed rms')
fdsplotlib.plot_to_fig(x_data=0.9*uprof, y_data=z, figure_handle=fig, marker_style='k--', data_label=None)
fdsplotlib.plot_to_fig(x_data=umean, y_data=z, figure_handle=fig, marker_style='b>-', data_label='FDS mean')

k = M.columns.get_loc('urms')
z = M.iloc[:, k-1].values
urms = M.iloc[:, k].values

fdsplotlib.plot_to_fig(x_data=umean+urms, y_data=z, figure_handle=fig, marker_style='r--', data_label='FDS rms')
fdsplotlib.plot_to_fig(x_data=umean-urms, y_data=z, figure_handle=fig, marker_style='r--', data_label=None)

plt.savefig(pltdir + 'sem_par_leddy_p2.pdf', format='pdf')
plt.close()

umean_error = np.linalg.norm(umean - uprof) / umax / len(umean)
if umean_error > error_tolerance:
    print(f'Matlab Warning: sem_par_leddy_p2.fds umean_error = {umean_error}')

urms_error = np.linalg.norm(umean + urms - 1.1*uprof) / umax / len(umean)
if urms_error > error_tolerance:
    print(f'Matlab Warning: sem_par_leddy_p2.fds urms_error = {urms_error}')

# Atmospheric profile

M = pd.read_csv(outdir + 'sem_atm_leddy_p2_line.csv', header=1)

k = M.columns.get_loc('umean')
z = M.iloc[:, k-1].values
umean = M.iloc[:, k].values

u0 = 1
z0 = 0.5
p = 0.3
uprof = u0 * (z / z0)**p

fig = fdsplotlib.plot_to_fig(x_data=uprof, y_data=z, marker_style='k-', data_label='Prescribed mean',
                             revision_label=version_string,
                             x_label='$u$ (m/s)', y_label='$z$ (m)',
                             x_min=0, x_max=1.4, y_min=0, y_max=1)

fdsplotlib.plot_to_fig(x_data=1.1*uprof, y_data=z, figure_handle=fig, marker_style='k--', data_label='Prescribed rms')
fdsplotlib.plot_to_fig(x_data=0.9*uprof, y_data=z, figure_handle=fig, marker_style='k--', data_label=None)
fdsplotlib.plot_to_fig(x_data=umean, y_data=z, figure_handle=fig, marker_style='b>-', data_label='FDS mean')

k = M.columns.get_loc('urms')
z = M.iloc[:, k-1].values
urms = M.iloc[:, k].values

fdsplotlib.plot_to_fig(x_data=umean+urms, y_data=z, figure_handle=fig, marker_style='r--', data_label='FDS rms')
fdsplotlib.plot_to_fig(x_data=umean-urms, y_data=z, figure_handle=fig, marker_style='r--', data_label=None)

plt.savefig(pltdir + 'sem_atm_leddy_p2.pdf', format='pdf')
plt.close()

umean_error = np.linalg.norm(umean - uprof) / u0 / len(umean)
if umean_error > error_tolerance:
    print(f'Matlab Warning: sem_atm_leddy_p2.fds umean_error = {umean_error}')

urms_error = np.linalg.norm(umean + urms - 1.1*uprof) / u0 / len(umean)
if urms_error > error_tolerance:
    print(f'Matlab Warning: sem_atm_leddy_p2.fds urms_error = {urms_error}')

# RAMP profile

M = pd.read_csv(outdir + 'sem_ramp_leddy_p2_line.csv', header=1)

k = M.columns.get_loc('umean')
z = M.iloc[:, k-1].values
umean = M.iloc[:, k].values

u0 = 1
uprof = np.zeros(len(z))
for i in range(len(z)):
    if z[i] < 0.5:
        uprof[i] = z[i] * 2 * u0
    else:
        uprof[i] = (1 - z[i]) * 2 * u0

fig = fdsplotlib.plot_to_fig(x_data=uprof, y_data=z, marker_style='k-', data_label='Prescribed mean',
                             revision_label=version_string,
                             x_label='$u$ (m/s)', y_label='$z$ (m)',
                             x_min=0, x_max=1.2, y_min=0, y_max=1)

fdsplotlib.plot_to_fig(x_data=1.1*uprof, y_data=z, figure_handle=fig, marker_style='k--', data_label='Prescribed rms')
fdsplotlib.plot_to_fig(x_data=0.9*uprof, y_data=z, figure_handle=fig, marker_style='k--', data_label=None)
fdsplotlib.plot_to_fig(x_data=umean, y_data=z, figure_handle=fig, marker_style='b>-', data_label='FDS mean')

k = M.columns.get_loc('urms')
z = M.iloc[:, k-1].values
urms = M.iloc[:, k].values

fdsplotlib.plot_to_fig(x_data=umean+urms, y_data=z, figure_handle=fig, marker_style='r--', data_label='FDS rms')
fdsplotlib.plot_to_fig(x_data=umean-urms, y_data=z, figure_handle=fig, marker_style='r--', data_label=None)

plt.savefig(pltdir + 'sem_ramp_leddy_p2.pdf', format='pdf')
plt.close()

umean_error = np.linalg.norm(umean - uprof) / u0 / len(umean)
if umean_error > error_tolerance:
    print(f'Matlab Warning: sem_ramp_leddy_p2.fds umean_error = {umean_error}')

urms_error = np.linalg.norm(umean + urms - 1.1*uprof) / u0 / len(umean)
if urms_error > error_tolerance:
    print(f'Matlab Warning: sem_ramp_leddy_p2.fds urms_error = {urms_error}')

# Monin-Obukhov profile at OPEN inflow boundary VELOCITY

# Expected values
E = pd.read_csv(outdir + 'sem_open_wind_MO_profile.csv', header=0)
z_exp = E['z'].values
u_exp = E['u'].values
T_exp = E['T'].values

M = pd.read_csv(outdir + 'sem_open_wind_line.csv', header=1)

z_fds = M['UMEAN-z'].values
u_fds = M['UMEAN'].values
u_fds_rms = M['URMS'].values
T_fds = M['TMEAN'].values
T_fds_rms = M['TRMS'].values

I = 0.1  # turbulence intensity (I), VEL_RMS=1 m/s, U_REF = 10 m/s from input file

fig = fdsplotlib.plot_to_fig(x_data=u_exp, y_data=z_exp, marker_style='k>', data_label='Monin-Obukhov profile',
                             revision_label=version_string,
                             x_label='$u$ (m/s)', y_label='$z$ (m)',
                             x_min=0, x_max=12, y_min=0, y_max=10)

fdsplotlib.plot_to_fig(x_data=(1+I)*u_exp, y_data=z_exp, figure_handle=fig, marker_style='k:', data_label='Prescribed rms')
fdsplotlib.plot_to_fig(x_data=(1-I)*u_exp, y_data=z_exp, figure_handle=fig, marker_style='k:', data_label=None)
fdsplotlib.plot_to_fig(x_data=u_fds, y_data=z_fds, figure_handle=fig, marker_style='b-', data_label='FDS mean')
fdsplotlib.plot_to_fig(x_data=u_fds+u_fds_rms, y_data=z_fds, figure_handle=fig, marker_style='b--', data_label='FDS rms')
fdsplotlib.plot_to_fig(x_data=u_fds-u_fds_rms, y_data=z_fds, figure_handle=fig, marker_style='b--', data_label=None)

plt.savefig(pltdir + 'sem_open_wind_u_prof.pdf', format='pdf')
plt.close()

u_ref = 10
kk = np.where((z_exp < np.max(z_fds)) & (z_exp > np.min(z_fds)))[0]
u_fds_int = np.interp(z_exp[kk], z_fds, u_fds)
umean_error = np.linalg.norm(u_exp[kk] - u_fds_int) / u_ref / len(u_fds_int)
if umean_error > error_tolerance:
    print(f'Matlab Warning: sem_open_wind.fds umean_error = {umean_error}')

u_fds_rms_int = np.interp(z_exp[kk], z_fds, u_fds_rms)
urms_error = np.linalg.norm(u_fds_int + u_fds_rms_int - (1+I)*u_exp[kk]) / u_ref / len(u_fds_rms_int)
if urms_error > error_tolerance:
    print(f'Matlab Warning: sem_open_wind.fds urms_error = {urms_error}')

# Monin-Obukhov profile at OPEN inflow boundary TEMPERATURE

fig = fdsplotlib.plot_to_fig(x_data=T_exp, y_data=z_exp, marker_style='ko', data_label='Monin-Obukhov profile',
                             revision_label=version_string,
                             x_min=19, x_max=20, y_min=0, y_max=10,
                             x_label=r'$T$ ($^\circ$C)', y_label='$z$ (m)')

fdsplotlib.plot_to_fig(x_data=T_fds, y_data=z_fds, figure_handle=fig, marker_style='r-', data_label='FDS mean')

plt.savefig(pltdir + 'sem_open_wind_T_prof.pdf', format='pdf')
plt.close()

T_ref = 20 + 273
kk = np.where((z_exp < np.max(z_fds)) & (z_exp > np.min(z_fds)))[0]
T_fds_int = np.interp(z_exp[kk], z_fds, T_fds)
Tmean_error = np.linalg.norm(T_exp[kk] - T_fds_int) / T_ref / len(T_fds_int)
if Tmean_error > error_tolerance:
    print(f'Matlab Warning: sem_open_wind.fds Tmean_error = {Tmean_error}')

