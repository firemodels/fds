
# FDS Verification Guide, Moody Chart

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
import os
import fdsplotlib

plot_style = fdsplotlib.get_plot_style('fds')

outdir = '../../../out/Moody_Chart/'
pltdir = '../../Manuals/FDS_Verification_Guide/SCRIPT_FIGURES/'

git_file = outdir + 'poiseuille_N64_mu025_git.txt'
version_string = fdsplotlib.get_version_string(git_file)

def colebrook(Re, RR, ff, tol):
    """
    Solve Colebrook equation for friction factor

    Parameters:
    Re - Reynolds number
    RR - Relative roughness (epsilon/D)
    ff - Initial guess for friction factor
    tol - Tolerance for convergence

    Returns:
    f_cor - Corrected friction factor
    error - Final error
    iter - Number of iterations
    """
    max_iter = 1000
    f = ff
    error = 1.0
    iter = 0

    # Colebrook-White equation: 1/sqrt(f) = -2*log10(RR/3.7 + 2.51/(Re*sqrt(f)))
    while error > tol and iter < max_iter:
        f_old = f
        if Re > 0 and f > 0:
            sqrt_f = np.sqrt(f)
            term1 = RR / 3.7
            term2 = 2.51 / (Re * sqrt_f)
            f_new = 1.0 / (-2.0 * np.log10(term1 + term2))**2
            f = f_new
            error = abs(f - f_old)
        else:
            break
        iter += 1

    f_cor = f
    return f_cor, error, iter


def friction_factor_calc(dpdx, H, filename, *args):
    """
    Calculate friction factor and Reynolds number from CSV data

    Parameters:
    -----------
    dpdx : float
        Pressure gradient
    H : float
        Channel height
    filename : str
        CSV filename to read
    *args : optional
        Optional viscosity value (mu)

    Returns:
    --------
    f_fds : float
        Friction factor from FDS
    Re_H : float
        Reynolds number based on H
    """

    # Read CSV file, skipping first 2 rows
    M = np.genfromtxt(filename, delimiter=',', skip_header=2)

    # Extract velocity data from column 2
    ubar = M[:, 1]

    # Handle optional viscosity argument
    if len(args) == 0:
        # If no optional argument provided (nargin==3), use max of column 4
        mu = np.max(M[:, 3])
    elif len(args) == 1:
        # If optional argument provided (nargin==4), use it
        mu = args[0]
    else:
        raise ValueError("Too many arguments provided")

    # Extract density from column 5 (Matlab 1-indexed column 5 = Python 0-indexed column 4)
    rho = np.max(M[:, 4])

    # Steady state mean velocity (planar averaged)
    U = ubar[-1]  # Last element of ubar array

    # Reynolds number based on H
    Re_H = H * U * rho / mu

    # Friction factor from FDS
    f_fds = 2 * (-dpdx) * H / (rho * U**2)

    return f_fds, Re_H


n = 100
Re = np.logspace(3.3, 8, n)
RR = [0, 1e-4, 1e-3, 1e-2, 1e-1]
tol = 1e-3

# Initialize figure
fig = fdsplotlib.plot_to_fig(x_data=[1e-3,1e-3], y_data=[1e-3,1e-3],
                             x_min=1e2, x_max=1e8, y_min=0.005, y_max=0.2,
                             plot_type='loglog',
                             plot_origin=(0.8,plot_style['Plot_Y']),
                             revision_label=version_string,
                             x_label='Re$_H$',
                             y_label='$f$')

# Initialize storage
M = np.zeros((n, len(RR) + 1))
f = np.zeros(n)

for i in range(len(RR)):

    for j in range(n):
        # initial guess for f
        if j > 0:
            ff = f[j-1]
        else:
            ff = 0.1

        f[j], error, iter = colebrook(Re[j], RR[i], ff, tol)
        M[j, 0] = Re[j]
        M[j, i+1] = f[j]

    if i == 0:
        fdsplotlib.plot_to_fig(x_data=Re, y_data=f, figure_handle=fig, marker_style='k-')
    if i == 1:
        fdsplotlib.plot_to_fig(x_data=Re, y_data=f, figure_handle=fig, marker_style='b-')
    if i == 2:
        fdsplotlib.plot_to_fig(x_data=Re, y_data=f, figure_handle=fig, marker_style='r-')
    if i == 3:
        fdsplotlib.plot_to_fig(x_data=Re, y_data=f, figure_handle=fig, marker_style='g-')
    if i == 4:
        fdsplotlib.plot_to_fig(x_data=Re, y_data=f, figure_handle=fig, marker_style='m-')

Re_DNS = np.logspace(2, 3.3)
f_DNS = 24.0 / Re_DNS
fdsplotlib.plot_to_fig(x_data=Re_DNS, y_data=f_DNS, figure_handle=fig, marker_style='k-')

ax = plt.gca()
ax.set_xlabel(r'Re$_H$', fontsize=plot_style['Label_Font_Size'])
ax.set_ylabel(r'$f$', fontsize=plot_style['Label_Font_Size'], rotation=0.0)

# gather FDS results (laminar)
L = 1
dpdx = -1
f_fds = np.zeros(2)
Re_fds = np.zeros(2)
f_fds[0], Re_fds[0] = friction_factor_calc(dpdx, L, outdir + 'poiseuille_N64_mu025_devc.csv')
f_fds[1], Re_fds[1] = friction_factor_calc(dpdx, L, outdir + 'poiseuille_N64_mu0125_devc.csv')

fdsplotlib.plot_to_fig(x_data=Re_fds, y_data=f_fds, figure_handle=fig, marker_style='k*')

# gather FDS results (turbulent)
mu = 1.84e-5

dpdx = -0.01
f, Re = friction_factor_calc(dpdx, L, outdir + 'moody_dpdx=-0p01_N8_devc.csv', mu)
fdsplotlib.plot_to_fig(x_data=Re, y_data=f, figure_handle=fig, marker_style='ksq')
f, Re = friction_factor_calc(dpdx, L, outdir + 'moody_dpdx=-0p01_N16_devc.csv', mu)
fdsplotlib.plot_to_fig(x_data=Re, y_data=f, figure_handle=fig, marker_style='k^')
f, Re = friction_factor_calc(dpdx, L, outdir + 'moody_dpdx=-0p01_N32_devc.csv', mu)
fdsplotlib.plot_to_fig(x_data=Re, y_data=f, figure_handle=fig, marker_style='ko')

dpdx = -1
f, Re = friction_factor_calc(dpdx, L, outdir + 'moody_dpdx=-1_N8_devc.csv', mu)
fdsplotlib.plot_to_fig(x_data=Re, y_data=f, figure_handle=fig, marker_style='ksq')
f, Re = friction_factor_calc(dpdx, L, outdir + 'moody_dpdx=-1_N16_devc.csv', mu)
fdsplotlib.plot_to_fig(x_data=Re, y_data=f, figure_handle=fig, marker_style='k^')
f, Re = friction_factor_calc(dpdx, L, outdir + 'moody_dpdx=-1_N32_devc.csv', mu)
fdsplotlib.plot_to_fig(x_data=Re, y_data=f, figure_handle=fig, marker_style='ko')

dpdx = -100
f, Re = friction_factor_calc(dpdx, L, outdir + 'moody_dpdx=-100_N8_devc.csv', mu)
fdsplotlib.plot_to_fig(x_data=Re, y_data=f, figure_handle=fig, marker_style='ksq', data_label='$N_z=8$, $H=1$')
f, Re = friction_factor_calc(dpdx, L, outdir + 'moody_dpdx=-100_N16_devc.csv', mu)
fdsplotlib.plot_to_fig(x_data=Re, y_data=f, figure_handle=fig, marker_style='k^', data_label='$N_z=16$, $H=1$')
f, Re = friction_factor_calc(dpdx, L, outdir + 'moody_dpdx=-100_N32_devc.csv', mu)
fdsplotlib.plot_to_fig(x_data=Re, y_data=f, figure_handle=fig, marker_style='ko', data_label='$N_z=32$, $H=1$')

dpdx = -0.0001
f, Re = friction_factor_calc(dpdx, L, outdir + 's=p0001_dpdx=-p0001_N8_devc.csv', mu)
fdsplotlib.plot_to_fig(x_data=Re, y_data=f, figure_handle=fig, marker_style='bsq')
f, Re = friction_factor_calc(dpdx, L, outdir + 's=p001_dpdx=-p0001_N8_devc.csv', mu)
fdsplotlib.plot_to_fig(x_data=Re, y_data=f, figure_handle=fig, marker_style='rsq')
f, Re = friction_factor_calc(dpdx, L, outdir + 's=p01_dpdx=-p0001_N8_devc.csv', mu)
fdsplotlib.plot_to_fig(x_data=Re, y_data=f, figure_handle=fig, marker_style='gsq')
f, Re = friction_factor_calc(dpdx, L, outdir + 's=p1_dpdx=-p0001_N8_devc.csv', mu)
fdsplotlib.plot_to_fig(x_data=Re, y_data=f, figure_handle=fig, marker_style='msq')

dpdx = -0.01
f, Re = friction_factor_calc(dpdx, L, outdir + 's=p0001_dpdx=-p01_N8_devc.csv', mu)
fdsplotlib.plot_to_fig(x_data=Re, y_data=f, figure_handle=fig, marker_style='bsq')
f, Re = friction_factor_calc(dpdx, L, outdir + 's=p001_dpdx=-p01_N8_devc.csv', mu)
fdsplotlib.plot_to_fig(x_data=Re, y_data=f, figure_handle=fig, marker_style='rsq')
f, Re = friction_factor_calc(dpdx, L, outdir + 's=p01_dpdx=-p01_N8_devc.csv', mu)
fdsplotlib.plot_to_fig(x_data=Re, y_data=f, figure_handle=fig, marker_style='gsq')
f, Re = friction_factor_calc(dpdx, L, outdir + 's=p1_dpdx=-p01_N8_devc.csv', mu)
fdsplotlib.plot_to_fig(x_data=Re, y_data=f, figure_handle=fig, marker_style='msq')

dpdx = -1
f, Re = friction_factor_calc(dpdx, L, outdir + 's=p0001_dpdx=-1_N8_devc.csv', mu)
fdsplotlib.plot_to_fig(x_data=Re, y_data=f, figure_handle=fig, marker_style='bsq')
f, Re = friction_factor_calc(dpdx, L, outdir + 's=p001_dpdx=-1_N8_devc.csv', mu)
fdsplotlib.plot_to_fig(x_data=Re, y_data=f, figure_handle=fig, marker_style='rsq')
f, Re = friction_factor_calc(dpdx, L, outdir + 's=p01_dpdx=-1_N8_devc.csv', mu)
fdsplotlib.plot_to_fig(x_data=Re, y_data=f, figure_handle=fig, marker_style='gsq')
f, Re = friction_factor_calc(dpdx, L, outdir + 's=p1_dpdx=-1_N8_devc.csv', mu)
fdsplotlib.plot_to_fig(x_data=Re, y_data=f, figure_handle=fig, marker_style='msq')
f, Re = friction_factor_calc(dpdx, L, outdir + 's=p0001_dpdx=-1_N16_devc.csv', mu)
fdsplotlib.plot_to_fig(x_data=Re, y_data=f, figure_handle=fig, marker_style='b^')
f, Re = friction_factor_calc(dpdx, L, outdir + 's=p001_dpdx=-1_N16_devc.csv', mu)
fdsplotlib.plot_to_fig(x_data=Re, y_data=f, figure_handle=fig, marker_style='r^')
f, Re = friction_factor_calc(dpdx, L, outdir + 's=p01_dpdx=-1_N16_devc.csv', mu)
fdsplotlib.plot_to_fig(x_data=Re, y_data=f, figure_handle=fig, marker_style='g^')
f, Re = friction_factor_calc(dpdx, 2*L, outdir + 's=p02_dpdx=-1_N16_devc.csv', mu)
fdsplotlib.plot_to_fig(x_data=Re, y_data=f, figure_handle=fig, marker_style='g>', data_label='$N_z=16$, $H=2$')

dpdx = -100
f, Re = friction_factor_calc(dpdx, L, outdir + 's=p0001_dpdx=-100_N8_devc.csv', mu)
fdsplotlib.plot_to_fig(x_data=Re, y_data=f, figure_handle=fig, marker_style='bsq')
f, Re = friction_factor_calc(dpdx, L, outdir + 's=p001_dpdx=-100_N8_devc.csv', mu)
fdsplotlib.plot_to_fig(x_data=Re, y_data=f, figure_handle=fig, marker_style='rsq')
f, Re = friction_factor_calc(dpdx, L, outdir + 's=p01_dpdx=-100_N8_devc.csv', mu)
fdsplotlib.plot_to_fig(x_data=Re, y_data=f, figure_handle=fig, marker_style='gsq')
f, Re = friction_factor_calc(dpdx, L, outdir + 's=p1_dpdx=-100_N8_devc.csv', mu)
fdsplotlib.plot_to_fig(x_data=Re, y_data=f, figure_handle=fig, marker_style='msq')
f, Re = friction_factor_calc(dpdx, L, outdir + 's=p0001_dpdx=-100_N16_devc.csv', mu)
fdsplotlib.plot_to_fig(x_data=Re, y_data=f, figure_handle=fig, marker_style='b^')
f, Re = friction_factor_calc(dpdx, L, outdir + 's=p001_dpdx=-100_N16_devc.csv', mu)
fdsplotlib.plot_to_fig(x_data=Re, y_data=f, figure_handle=fig, marker_style='r^')
f, Re = friction_factor_calc(dpdx, L, outdir + 's=p01_dpdx=-100_N16_devc.csv', mu)
fdsplotlib.plot_to_fig(x_data=Re, y_data=f, figure_handle=fig, marker_style='g^')
f, Re = friction_factor_calc(dpdx, 2*L, outdir + 's=p02_dpdx=-100_N16_devc.csv', mu)
fdsplotlib.plot_to_fig(x_data=Re, y_data=f, figure_handle=fig, marker_style='g>')

plt.text(1.3e8, 2e-1, r'$s/H$', fontsize=plot_style['Key_Font_Size'])
plt.text(1.3e8, 1.2e-2, str(RR[1]), fontsize=plot_style['Key_Font_Size'])
plt.text(1.3e8, 1.95e-2, str(RR[2]), fontsize=plot_style['Key_Font_Size'])
plt.text(1.3e8, 3.85e-2, str(RR[3]), fontsize=plot_style['Key_Font_Size'])
plt.text(1.3e8, 1.02e-1, str(RR[4]), fontsize=plot_style['Key_Font_Size'])

plt.savefig(pltdir + 'fds_moody_chart.pdf', format='pdf')
plt.close()


# Poiseuille convergence plot

dpdx = -1
L = 1
N = [8, 16, 32, 64]

f = [0.0] * 4
Re = [0.0] * 4

f[0], Re[0] = friction_factor_calc(dpdx, L, outdir + 'poiseuille_N8_mu025_devc.csv')
f[1], Re[1] = friction_factor_calc(dpdx, L, outdir + 'poiseuille_N16_mu025_devc.csv')
f[2], Re[2] = friction_factor_calc(dpdx, L, outdir + 'poiseuille_N32_mu025_devc.csv')
f[3], Re[3] = friction_factor_calc(dpdx, L, outdir + 'poiseuille_N64_mu025_devc.csv')

f = np.array(f)
Re = np.array(Re)
N = np.array(N)

dz = L / N
error = np.abs(f - 24.0 / Re)

fig = fdsplotlib.plot_to_fig(x_data=dz, y_data=error, marker_style='b*-', data_label='FDS',
                             x_min=0.01, x_max=0.2, y_min=5e-5, y_max=0.01,
                             plot_type='loglog',
                             revision_label=version_string,
                             x_label=r'Grid Spacing, $\delta z$ (m)',
                             y_label='Friction Factor Error')

fdsplotlib.plot_to_fig(x_data=dz, y_data=0.05*dz,   figure_handle=fig, marker_style='k--', data_label=r'$O(\delta z)$')
fdsplotlib.plot_to_fig(x_data=dz, y_data=0.4*dz**2, figure_handle=fig, marker_style='k-' , data_label=r'$O(\delta z^2)$')

plt.savefig(pltdir + 'poiseuille_convergence.pdf', format='pdf')
plt.close()


# Poiseuille convergence plot using complex geometry

dpdx = -1
L = 1
N = np.array([10, 20, 40, 80])

outdir = '../../Verification/Complex_Geometry/'

# Method for cross velocity forcing:
vmethod = ['_stm']

for im in range(1):
    
    mth = vmethod[im]
    
    # plot convergence for Poiseuille flow aligned case theta=0 (mu = 0.025)

    f = np.zeros(4)
    Re = np.zeros(4)
    f2 = np.zeros(4)
    Re2 = np.zeros(4)
    
    f[0], Re[0] = friction_factor_calc(dpdx, L, outdir + 'geom_poiseuille_N10a_theta0' + mth + '_devc.csv')
    f[1], Re[1] = friction_factor_calc(dpdx, L, outdir + 'geom_poiseuille_N20a_theta0' + mth + '_devc.csv')
    f[2], Re[2] = friction_factor_calc(dpdx, L, outdir + 'geom_poiseuille_N40a_theta0' + mth + '_devc.csv')
    f[3], Re[3] = friction_factor_calc(dpdx, L, outdir + 'geom_poiseuille_N80a_theta0' + mth + '_devc.csv')
    
    f2[0], Re2[0] = friction_factor_calc(dpdx, L, outdir + 'geom_poiseuille_N10nah_theta0' + mth + '_devc.csv')
    f2[1], Re2[1] = friction_factor_calc(dpdx, L, outdir + 'geom_poiseuille_N20nah_theta0' + mth + '_devc.csv')
    f2[2], Re2[2] = friction_factor_calc(dpdx, L, outdir + 'geom_poiseuille_N40nah_theta0' + mth + '_devc.csv')
    f2[3], Re2[3] = friction_factor_calc(dpdx, L, outdir + 'geom_poiseuille_N80nah_theta0' + mth + '_devc.csv')
    
    dz = L / N
    error = np.abs(f - 24.0 / Re)
    error2 = np.abs(f2 - 24.0 / Re2)
    
    fig = fdsplotlib.plot_to_fig(x_data=dz, y_data=error, marker_style='b*-', data_label='FDS, $h=0$',
                                 x_min=0.01, x_max=0.2, y_min=5e-5, y_max=0.01,
                                 plot_type='loglog',
                                 revision_label=version_string,
                                 x_label=r'Grid Spacing, $\delta z$ (m)',
                                 y_label='Friction Factor Error')

    fdsplotlib.plot_to_fig(x_data=dz, y_data=error2,   figure_handle=fig, marker_style='rx-', data_label=r'$h=\delta z/3$')
    fdsplotlib.plot_to_fig(x_data=dz, y_data=0.12*dz,  figure_handle=fig, marker_style='k--', data_label=r'$O(\delta z)$')
    fdsplotlib.plot_to_fig(x_data=dz, y_data=0.4*dz**2,figure_handle=fig, marker_style='k-' , data_label=r'$O(\delta z^2)$')
    
    output_file = pltdir + 'geom_poiseuille_convergence_theta0a' + mth + '.pdf'
    plt.savefig(output_file, format='pdf')
    plt.close()
    
    # plot convergence for Poiseuille flow not aligned case theta=0 (mu = 0.025)
    
    f = np.zeros(4)
    Re = np.zeros(4)
    f2 = None
    Re2 = None
    H = []
    
    f[0], Re[0] = friction_factor_calc(dpdx, L, outdir + 'geom_poiseuille_N10na_theta0' + mth + '_devc.csv')
    f[1], Re[1] = friction_factor_calc(dpdx, L, outdir + 'geom_poiseuille_N20na_theta0' + mth + '_devc.csv')
    f[2], Re[2] = friction_factor_calc(dpdx, L, outdir + 'geom_poiseuille_N40na_theta0' + mth + '_devc.csv')
    f[3], Re[3] = friction_factor_calc(dpdx, L, outdir + 'geom_poiseuille_N80na_theta0' + mth + '_devc.csv')
    
    dz = L / N
    error = np.abs(f - 24.0 / Re)
    
    fig = fdsplotlib.plot_to_fig(x_data=dz, y_data=error, marker_style='b*-', data_label=r'FDS, $h=\delta z_{10}/11$',
                                 x_min=0.01, x_max=0.2, y_min=5e-5, y_max=0.01,
                                 plot_type='loglog',
                                 revision_label=version_string,
                                 x_label=r'Grid Spacing, $\delta z$ (m)',
                                 y_label='Friction Factor Error')

    fdsplotlib.plot_to_fig(x_data=dz, y_data=0.05*dz,  figure_handle=fig, marker_style='k--', data_label=r'$O(\delta z)$')
    fdsplotlib.plot_to_fig(x_data=dz, y_data=0.4*dz**2,figure_handle=fig, marker_style='k-' , data_label=r'$O(\delta z^2)$')
    
    output_file = pltdir + 'geom_poiseuille_convergence_theta0na' + mth + '.pdf'
    plt.savefig(output_file, format='pdf')
    plt.close()


