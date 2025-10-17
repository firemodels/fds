
# Heat transfer within a sphere
# The solution is from Carslaw and Jaeger, Sec. 9.8, p. 243, Eq. (6)
# The notation adopted here follows C. Lautenberger, IAFSS, 2014.

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
import sys
import fdsplotlib 

plot_style = fdsplotlib.get_plot_style('fds')

outdir = '../../Verification/Heat_Transfer/'
pltdir = '../../Manuals/FDS_Verification_Guide/SCRIPT_FIGURES/'

git_file = outdir + 'ht3d_sphere_48_git.txt'
version_string = fdsplotlib.get_version_string(git_file)


# analytical solution

k = 1.0          # W/m/k
rho = 1000       # kg/m3
cp = 1000        # J/kg/K
g0 = 2e5         # W/m3
alpha = k/(rho*cp)  # m2/s

a1 = 0.1  # m (this should match the last cell face in ht3d_sphere_96.fds)
a2 = 0.1  # m (this should match the last cell face in ht3d_sphere_48.fds)
a3 = 0.1  # m (this should match the last cell face in ht3d_sphere_24.fds)

n1 = 41
n2 = 21
n3 = 11
r1 = np.linspace(0.00001, 0.1, n1)  # this should match line DEVC in ht3d_sphere_96.fds
r2 = np.linspace(0.00001, 0.1, n2)  # this should match line DEVC in ht3d_sphere_48.fds
r3 = np.linspace(0.00001, 0.1, n3)  # this should match line DEVC in ht3d_sphere_24.fds

t = [10, 20, 60, 120, 180]  # seconds
markers = ['bo','go','ro','co','mo']

fig = fdsplotlib.plot_to_fig(x_data=[-1,-1], y_data=[-1,-1],
                             x_min=0, x_max=0.105, y_min=20, y_max=60,
                             revision_label=version_string,
                             x_label='Radial Distance (m)',
                             y_label=r'Temperature ($^\circ$C)')

for m in range(len(t)):

    DT1 = np.zeros(len(r1))
    for ii in range(len(r1)):
        sum_term = 0
        for n in range(1, n1+1):
            sum_term = sum_term + (-1)**n/n**3 * np.sin(n*np.pi*r1[ii]/a1) * np.exp(-alpha*t[m]*(n*np.pi/a1)**2)
        DT1[ii] = 20 + g0/(6*k) * (a1**2 - r1[ii]**2) + 2*g0*a1**3/(k*np.pi**3*r1[ii]) * sum_term

    DT2 = np.zeros(len(r2))
    for jj in range(len(r2)):
        sum_term = 0
        for n in range(1, n2+1):
            sum_term = sum_term + (-1)**n/n**3 * np.sin(n*np.pi*r2[jj]/a2) * np.exp(-alpha*t[m]*(n*np.pi/a2)**2)
        DT2[jj] = 20 + g0/(6*k) * (a2**2 - r2[jj]**2) + 2*g0*a2**3/(k*np.pi**3*r2[jj]) * sum_term

    if m==0:
        fdsplotlib.plot_to_fig(x_data=r2, y_data=DT2, figure_handle=fig, marker_style=markers[m], data_label='Exact')
    else:
        fdsplotlib.plot_to_fig(x_data=r2, y_data=DT2, figure_handle=fig, marker_style=markers[m])

    DT3 = np.zeros(len(r3))
    for kk in range(len(r3)):
        sum_term = 0
        for n in range(1, n3+1):
            sum_term = sum_term + (-1)**n/n**3 * np.sin(n*np.pi*r3[kk]/a3) * np.exp(-alpha*t[m]*(n*np.pi/a3)**2)
        DT3[kk] = 20 + g0/(6*k) * (a3**2 - r3[kk]**2) + 2*g0*a3**3/(k*np.pi**3*r3[kk]) * sum_term

# gather FDS results

fnt = ['ht3d_sphere_48']
fileName = ['ht3d_sphere_24', 'ht3d_sphere_48', 'ht3d_sphere_96']
nc_array = np.array([25, 50, 100])
dx_array = 0.25 / nc_array

M = pd.read_csv(outdir + fnt[0] + '_prof_1.csv', sep=',', skiprows=3, header=None)
T_fds1 = M.iloc[ 1, 63:].values 
T_fds2 = M.iloc[ 2, 63:].values 
T_fds3 = M.iloc[ 6, 63:].values 
T_fds4 = M.iloc[12, 63:].values 
T_fds5 = M.iloc[-1, 63:].values
fdsplotlib.plot_to_fig(x_data=r2, y_data=T_fds1, figure_handle=fig, marker_style='b--', data_label='FDS $t=10$ s')
fdsplotlib.plot_to_fig(x_data=r2, y_data=T_fds2, figure_handle=fig, marker_style='g--', data_label='FDS $t=20$ s')
fdsplotlib.plot_to_fig(x_data=r2, y_data=T_fds3, figure_handle=fig, marker_style='r--', data_label='FDS $t=60$ s')
fdsplotlib.plot_to_fig(x_data=r2, y_data=T_fds4, figure_handle=fig, marker_style='c--', data_label='FDS $t=120$ s')
fdsplotlib.plot_to_fig(x_data=r2, y_data=T_fds5, figure_handle=fig, marker_style='m--', data_label='FDS $t=180$ s')

ERROR = abs(M.iloc[-1, 63] - DT2[1]) / (DT2[1] - 20)
if ERROR > 0.01:
    print(f'ERROR: ht3d_sphere cases out of tolerance. ERROR = {ERROR}')

plt.savefig(pltdir + 'ht3d_sphere_profile.pdf', format='pdf')
plt.close()

# estimating L1 & L2 norm errors

Linfe = []  # initialize L1 norm error vector
dxx = []    # init dxx vector
DT = [DT3[-1], DT2[-1], DT1[-1]]
r = [r3[-1], r2[-1], r1[-1]]

for i in range(len(fileName)):
    M1 = pd.read_csv(outdir + fileName[i] + '_prof_1.csv', sep=',', skiprows=3, header=0)
    T_fds = M1.iloc[-1, -1]  # FDS devices data (last row, last column)
    T_fds1 = T_fds
    Linfe.append(1/(1) * np.max(np.abs(DT[i] - T_fds1)))  # populates Linf norm error vector, element-by-element
    dxx.append(dx_array[i])  # populates dxx vector, element-by-element

Linfe = np.array(Linfe)
dxx = np.array(dxx)

fig = fdsplotlib.plot_to_fig(x_data=dxx, y_data=Linfe, marker_style='msq-', data_label='FDS',
                             x_min=2e-3, x_max=1e-2, y_min=5e-2, y_max=2,
                             plot_type='loglog',
                             revision_label=version_string,
                             x_label=r'$\delta x$ (m)',
                             y_label=r'$L_\infty$ Error ($^\circ$C)')

fdsplotlib.plot_to_fig(x_data=dxx, y_data=1e2*dxx,    marker_style='k--', data_label=r'$O(\delta x)$', figure_handle=fig)
fdsplotlib.plot_to_fig(x_data=dxx, y_data=1e4*dxx**2, marker_style='k-',  data_label=r'$O(\delta x^2)$', figure_handle=fig)

plt.savefig(pltdir + 'ht3d_sphere_convergence1.pdf', format='pdf')
plt.close()

