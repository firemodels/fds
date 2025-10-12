
# Plot comparisons of FDS vs Blasius boundary layer and Pohlausen thermal boundary layer

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import fdsplotlib

plot_style = fdsplotlib.get_plot_style('fds')

outdir = '../../Verification/Flowfields/'
pltdir = '../../Manuals/FDS_Verification_Guide/SCRIPT_FIGURES/'

git_file = outdir + 'blasius_16_git.txt'
version_string = fdsplotlib.get_version_string(git_file)


def blasius_analytic(u0, zmax, mu, rho, x, steps, fpp_init):
    """
    Parameters:
    u0 - free stream velocity
    zmax - maximum z coordinate for calculation
    mu - dynamic viscosity
    rho - density
    x - streamwise coordinate
    steps - number of steps between 0 and etamax
    fpp_init - initial (wall) value of 2nd derivative of Blasius function

    Returns:
    eta - the similarity coordinate normal to the wall
    f, fp - the Blasius function and its first derivative
    """

    etamax = zmax / np.sqrt(mu / rho * x / u0)
    deta = etamax / (steps - 1)
    eta = np.zeros(steps)
    f = np.zeros(steps)
    fp = np.zeros(steps)
    fpp = np.zeros(steps)
    fppp = np.zeros(steps)

    # initial guess for fpp
    fpp[0] = fpp_init

    for i in range(steps - 1):
        eta[i + 1] = eta[i] + deta

        # predictor
        # 1st
        k1 = np.zeros(3)
        k1[0] = fp[i]
        k1[1] = fpp[i]
        k1[2] = -f[i] * fpp[i] / 2

        fbar = f[i] + 0.5 * deta * k1[0]
        fpbar = fp[i] + 0.5 * deta * k1[1]
        fppbar = fpp[i] + 0.5 * deta * k1[2]

        # 2nd
        k2 = np.zeros(3)
        k2[0] = fpbar
        k2[1] = fppbar
        k2[2] = -fbar * fppbar / 2

        fbar = f[i] + 0.5 * deta * k2[0]
        fpbar = fp[i] + 0.5 * deta * k2[1]
        fppbar = fpp[i] + 0.5 * deta * k2[2]

        # 3rd
        k3 = np.zeros(3)
        k3[0] = fpbar
        k3[1] = fppbar
        k3[2] = -fbar * fppbar / 2

        fbar = f[i] + deta * k3[0]
        fpbar = fp[i] + deta * k3[1]
        fppbar = fpp[i] + deta * k3[2]

        # 4th
        k4 = np.zeros(3)
        k4[0] = fpbar
        k4[1] = fppbar
        k4[2] = -fbar * fppbar / 2

        # corrector
        f[i + 1] = f[i] + deta * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0]) / 6
        fp[i + 1] = fp[i] + deta * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1]) / 6
        fpp[i + 1] = fpp[i] + deta * (k1[2] + 2 * k2[2] + 2 * k3[2] + k4[2]) / 6
        fppp[i + 1] = -f[i + 1] * fpp[i + 1] / 2

    return eta, f, fp


M = pd.read_csv(outdir + 'blasius_16_line.csv', skiprows=2, header=None).values
u_16 = M[:, 1]
z_16 = M[:, 0]

M = pd.read_csv(outdir + 'blasius_32_line.csv', skiprows=2, header=None).values
u_32 = M[:, 1]
z_32 = M[:, 0]

M = pd.read_csv(outdir + 'blasius_64_line.csv', skiprows=2, header=None).values
u_64 = M[:, 1]
z_64 = M[:, 0]

# Plot Blasius profiles of FDS solutions

u0 = np.max(u_64)
zmax = 0.3
mu = 0.001
rho = 1.2
xc = 0.05

eta, f, fp = blasius_analytic(u0, zmax, mu, rho, xc, 257, 0.3318)
z_blasius = eta * np.sqrt(mu / rho * xc / u0)
u_blasius = fp * u0

range_indices = np.arange(0, len(u_blasius), 4)

fig = fdsplotlib.plot_to_fig(x_data=u_blasius[range_indices], y_data=z_blasius[range_indices], marker_style='ko', data_label='Blasius',
                             x_min=0, x_max=1.1, y_min=0, y_max=0.15,
                             revision_label=version_string,
                             x_label='$u$ (m/s)',
                             y_label='$z$ (m)')

fdsplotlib.plot_to_fig(x_data=u_16, y_data=z_16, figure_handle=fig, marker_style='k-.', data_label='$N_z=16$')
fdsplotlib.plot_to_fig(x_data=u_32, y_data=z_32, figure_handle=fig, marker_style='k--', data_label='$N_z=32$')
fdsplotlib.plot_to_fig(x_data=u_64, y_data=z_64, figure_handle=fig, marker_style='k-' , data_label='$N_z=64$')

plt.savefig(pltdir + 'blasius_profile.pdf', format='pdf')
plt.close()

# Error between FDS with Blasius solution

err = np.zeros(3)

#print(u_blasius)
for i in range(len(u_16)):
    err[0] = err[0] + (np.abs(u_16[i] - u_blasius[8 + i*16]))**2
err[0] = err[0] / 16
err[0] = np.sqrt(err[0])

for i in range(len(u_32)):
    err[1] = err[1] + (np.abs(u_32[i] - u_blasius[4 + i*8]))**2
err[1] = err[1] / 32
err[1] = np.sqrt(err[1])

for i in range(len(u_64)):
    err[2] = err[2] + (np.abs(u_64[i] - u_blasius[2 + i*4]))**2
err[2] = err[2] / 64
err[2] = np.sqrt(err[2])

dz = np.zeros(3)
dz[0] = np.abs(z_16[9] - z_16[8])
dz[1] = np.abs(z_32[9] - z_32[8])
dz[2] = np.abs(z_64[9] - z_64[8])

fig = fdsplotlib.plot_to_fig(x_data=dz, y_data=err, marker_style='k*-', data_label='FDS',
                             x_min=1e-3, x_max=1e-1, y_min=1e-3, y_max=1,
                             plot_type='loglog',
                             revision_label=version_string,
                             x_label=r'Grid Spacing, $\delta z$ (m)',
                             y_label='RMS Error (m/s)')

fdsplotlib.plot_to_fig(x_data=dz, y_data=10*dz    , figure_handle=fig, marker_style='k--', data_label=r'$O(\delta z)$')
fdsplotlib.plot_to_fig(x_data=dz, y_data=100*dz**2, figure_handle=fig, marker_style='k-',  data_label=r'$O(\delta z^2)$')

plt.savefig(pltdir + 'blasius_convergence.pdf', format='pdf')
plt.close()


# Pohlhausen solution, 1921.

outdir = '../../../out/Convection/'
chid = ['Pohlhausen_Pr_p5', 'Pohlhausen_Pr_1', 'Pohlhausen_Pr_2']
git_file = outdir + 'Pohlhausen_Pr_2_git.txt'
version_string = fdsplotlib.get_version_string(git_file)

# Blasius profile

Pr = [0.5, 1, 2]  # Prandtl numbers
u0 = 1.0          # see FDS input file
zmax = 1.0        # see FDS input file
mu = 0.001        # see FDS input and output files
rho = 1.2         # taken from FDS output file
nu = mu / rho
x = 5.0  # corresponds to the position of the DEVC array in the FDS input file
neta = 500
eta, f, fp = blasius_analytic(u0, zmax, mu, rho, x, neta, 0.3318)
etamax = eta[-1]

# thermal profile from Pohlhausen solution, Holman eq. (B-14)

fds_marker = ['k-','r-','g-']

fig = fdsplotlib.plot_to_fig(x_data=[-1,-1], y_data=[-1,-1], 
                             x_min=0, x_max=8, y_min=0, y_max=1,
                             revision_label=version_string,
                             x_label=r'$\eta=z/\sqrt{\nu x/u_0}$',
                             y_label=r'$\frac{T-T_0}{T_w-T_0}=1-\theta$')

for p in range(len(Pr)):

    deta = eta[1] - eta[0]

    I = np.where(eta <= etamax)[0]
    theta = np.zeros(len(I))

    for i in range(len(I)):
        ieta = I[i]
        num = 0.0
        for j in range(ieta + 1):
            intf = np.trapz(f[:j + 1]) * deta
            num = num + np.exp(-Pr[p] / 2.0 * intf) * deta
        denom = num
        for j in range(ieta + 1, neta):
            intf = np.trapz(f[:j + 1]) * deta
            denom = denom + np.exp(-Pr[p] / 2.0 * intf) * deta
        theta[i] = (num - deta) / (denom - deta)

    fdsplotlib.plot_to_fig(x_data=eta[I], y_data=1-theta, figure_handle=fig, marker_style=fds_marker[p], data_label=f'Pr={Pr[p]}')

plt.savefig(pltdir + 'Pohlhausen_similarity_solution.pdf', format='pdf')
plt.close()

# now plot the temperature profile in physical space

Tw = 21.0  # wall temperature
T0 = 20.0  # ambient temperature
fds_marker = ['ko-','ro-','go-']
p_marker   = ['k-','r-','g-']

fig = fdsplotlib.plot_to_fig(x_data=[-1,-1], y_data=[-1,-1], 
                             x_min=T0-0.1, x_max=Tw+0.1, y_min=0, y_max=1,
                             revision_label=version_string,
                             x_label=r'$T$ ($^\circ$C)',
                             y_label=r'$z$ (m)')

for p in range(len(Pr)):

    T = Tw * np.ones(len(theta)) + (T0 - Tw) * theta
    z = eta[I] * np.sqrt(nu * x / u0)

    fdsplotlib.plot_to_fig(x_data=T, y_data=z, figure_handle=fig, marker_style=p_marker[p], data_label=f'Pohlhausen Pr={Pr[p]}')

    M = pd.read_csv(outdir + chid[p] + '_line.csv', skiprows=1)
    zfds = M['z'].values
    Tfds = M['Tz'].values
    Ufds = M['Uz'].values

    fdsplotlib.plot_to_fig(x_data=Tfds[1:], y_data=zfds[1:], figure_handle=fig, marker_style=fds_marker[p], data_label=f'FDS Pr={Pr[p]}')

plt.savefig(pltdir + 'Pohlhausen_Tz_profile.pdf', format='pdf')
plt.close()

