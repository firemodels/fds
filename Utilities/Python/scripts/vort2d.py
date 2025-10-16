
# Generate all plots for the 2D Vortex section of the FDS Verification Guide

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import fdsplotlib

plot_style = fdsplotlib.get_plot_style('fds')

outdir = '../../Verification/NS_Analytical_Solution/'
pltdir = '../../Manuals/FDS_Verification_Guide/SCRIPT_FIGURES/'

git_file = outdir + 'vort2d_40_git.txt'
version_string = fdsplotlib.get_version_string(git_file)


L  = 0.3112
U0 = 35.0
Rc = 0.01556
BigGamma = 0.0359157

meshsize = [40, 80, 160, 320]
meshname = ['40', '80', '160', '320']
markers = ['k-','r-','c-','m-']
DT = 1.1114 * 10**(-4)
FLOWTIME = int(0.3112 / (U0 * DT))
TIMESTEPS = 4

DEVCNUM = np.zeros(4, dtype=int)
STEPSIZE = np.zeros(4)
COUNT = np.zeros((4, TIMESTEPS))
ERROR = np.zeros((4, TIMESTEPS))
M = {}

# Loop over four grid resolutions
for m in range(4): 

    DEVCNUM[m] = meshsize[m] // 5
    STEPSIZE[m] = 0.311200 / meshsize[m]

    M[m] = pd.read_csv(outdir + f'vort2d_{meshname[m]}_devc.csv', skiprows=2, header=None)

    # Plot u-velocity at x=0 Data for a Range of Times

    Z0 = STEPSIZE[m] / 2
    Z = np.arange(-2*Rc + Z0, 2*Rc - Z0 + STEPSIZE[m]/2, STEPSIZE[m])

    fig = fdsplotlib.plot_to_fig(x_data=U0-BigGamma*Z*np.exp(-Z**2/(2*Rc**2))/Rc**2, y_data=Z, marker_style='k-', data_label='Analytical',
                                 x_min=33.5, x_max=36.5, y_min=-0.03, y_max=0.03,
                                 revision_label=version_string,
                                 y_label='$z$ (m)', x_label='$u$ (m/s)')

    # Plot FDS simulation values with
    for t in range(TIMESTEPS):

        TIME = FLOWTIME * t
        COUNT[m, t] = 0.0
        ERROR[m, t] = 0.0
        UArray = np.zeros(DEVCNUM[m])

        for k in range(DEVCNUM[m]):
            UArray[k] = M[m].iloc[TIME, k+1]
            ZVAL = -2*Rc + Z0 + k*STEPSIZE[m]
            # Calculate analytical u-velocity values
            AxisExact = U0 - BigGamma * ZVAL * np.exp(-ZVAL**2/(2*Rc**2)) / Rc**2
            # Calculate initial error values
            ERROR[m, t] = ERROR[m, t] + ((UArray[k] - AxisExact))**2
            COUNT[m, t] = COUNT[m, t] + 1.0

        ERROR[m, t] = np.sqrt(ERROR[m, t] / COUNT[m, t])

        # Create plot of u-velocity along a line at x=0
        if t>0: fdsplotlib.plot_to_fig(x_data=UArray, y_data=Z, marker_style=markers[t], data_label=f'Pass {t}', figure_handle=fig)

    plt.savefig(pltdir + f'vort2d_{meshname[m]}_uzgraph.pdf', format='pdf')
    plt.close()

    # Plot u-Velocity at a point as a function of time

    T = np.arange(0.0, DT * (FLOWTIME * TIMESTEPS + 1) + DT, DT)
    XP = 0.0
    ZP = -2*Rc + Z0

    # Plot analytical solution (periodic to 6 loops)
    PointPlotExact = U0 - BigGamma * ZP * (
         np.exp(-(XP**2 + ZP**2 + (U0**2)*(T**2) -
                  2.*U0*XP*T) / (2.*Rc**2)) +
         np.exp(-((XP + 1*L)**2 + ZP**2 + (U0**2)*(T**2) -
                  2.*U0*(XP + 1*L)*T) / (2*Rc**2)) +
         np.exp(-((XP + 2*L)**2 + ZP**2 + (U0**2)*(T**2) -
                  2.*U0*(XP + 2*L)*T) / (2*Rc**2)) +
         np.exp(-((XP + 3*L)**2 + ZP**2 + (U0**2)*(T**2) -
                  2.*U0*(XP + 3*L)*T) / (2*Rc**2)) +
         np.exp(-((XP + 4*L)**2 + ZP**2 + (U0**2)*(T**2) -
                  2.*U0*(XP + 4*L)*T) / (2*Rc**2)) +
         np.exp(-((XP + 5*L)**2 + ZP**2 + (U0**2)*(T**2) -
                  2.*U0*(XP + 5*L)*T) / (2*Rc**2)) +
         np.exp(-((XP + 6*L)**2 + ZP**2 + (U0**2)*(T**2) -
                  2.*U0*(XP + 6*L)*T) / (2*Rc**2))) / Rc**2

    fig = fdsplotlib.plot_to_fig(x_data=T, y_data=PointPlotExact, marker_style='k-', data_label='Analytical',
                                 x_min=0, x_max=0.037, y_min=34.6, y_max=36,
                                 revision_label=version_string,
                                 x_label='Time (s)', y_label='$u$ (m/s)')

    T_fds = M[m].iloc[:, 0].values
    fdsplotlib.plot_to_fig(x_data=T_fds, y_data=M[m].iloc[:,1].values, figure_handle=fig, marker_style='r--', data_label='FDS')

    plt.savefig(pltdir + f'vort2d_{meshname[m]}_upgraph.pdf', format='pdf')
    plt.close()


# Generate Flow Diagram

U0   = 0.7
Min  = -0.05
Max  =  0.05
Step =  0.005

X, Z = np.meshgrid(np.arange(Min, Max + Step, Step), np.arange(Min, Max + Step, Step))
Psi = BigGamma * np.exp(-(X**2 + Z**2) / (2.0 * Rc**2))
DZ, DX = np.gradient(Psi, Step, Step)
U = U0 + DZ
W = -DX

fig = fdsplotlib.plot_to_fig(x_data=[-1,-1], y_data=[-1,-1], x_min=Min, x_max=Max, y_min=Min, y_max=Max,
                             figure_size=(plot_style['Scat_Paper_Width'],plot_style['Scat_Paper_Height']),
                             plot_size=(plot_style['Scat_Plot_Width'],plot_style['Scat_Plot_Height']),
                             plot_origin=(plot_style['Scat_Plot_X'],plot_style['Scat_Plot_Y']),
                             x_label='$x$ (m)', y_label='$z$ (m)',
                             revision_label=version_string)

ax = plt.gca()
ax.quiver(X, Z, U, W, scale=20, headwidth=5, headlength=7, headaxislength=6.5)

plt.savefig(pltdir + 'vort2d_diagram.pdf', format='pdf')
plt.close()


# Create convergence plot

dx = np.array([0.0077800, 0.0038900, 0.0019450, 0.0009725])
dx2 = np.array([0.0000605, 0.0000151, 0.0000038, 0.0000009])

fig = fdsplotlib.plot_to_fig(x_data=dx, y_data=120*dx, marker_style='k--', data_label=r'$O(\delta x)$',
                             x_min=0.0009, x_max=0.009, y_min=1e-3, y_max=1,
                             plot_type='loglog',
                             revision_label=version_string,
                             x_label='Grid Spacing (m)',
                             y_label='Error (m/s)')

fdsplotlib.plot_to_fig(x_data=dx, y_data=4500*dx2, marker_style='k-', data_label=r'$O(\delta x^2)$', figure_handle=fig)
fdsplotlib.plot_to_fig(x_data=dx, y_data=ERROR[0:4,1], marker_style='go-', data_label=r'FDS ($t=0.0089$ s)', figure_handle=fig)
fdsplotlib.plot_to_fig(x_data=dx, y_data=ERROR[0:4,2], marker_style='cv-', data_label=r'FDS ($t=0.0178$ s)', figure_handle=fig)
fdsplotlib.plot_to_fig(x_data=dx, y_data=ERROR[0:4,3], marker_style='r^-', data_label=r'FDS ($t=0.0267$ s)', figure_handle=fig)

plt.savefig(pltdir + 'vort2d_error.pdf', format='pdf')
plt.close()

