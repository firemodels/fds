
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import fdsplotlib

plot_style = fdsplotlib.get_plot_style('fds')

outdir = '../../Verification/Sprinklers_and_Sprays/'
pltdir = '../../Manuals/FDS_Verification_Guide/SCRIPT_FIGURES/'

git_file = outdir + 'terminal_velocity_dt_1_0_git.txt'
version_string = fdsplotlib.get_version_string(git_file)

M0 = pd.read_csv(outdir + 'terminal_velocity_dt_1_0_devc.csv', sep=',', header=1)
rhoa = M0.iloc[-1, 3]

g = 9.8
Cd = 1
rhod = 1000
D = 10e-3

tend = 20
eps = 1e-10
ttest = 0.1

K = 3 * rhoa * Cd / (4 * rhod * D)

dtvec = np.array([1.0, 0.1, 0.01, 0.001, 0.0001])
errvec = []
Linf = []
STIME = 20.
vtexact = np.sqrt(g / K)
zexact = lambda t: -np.log(np.cosh(np.sqrt(g * K) * t)) / K

M1 = pd.read_csv(outdir + 'terminal_velocity_dt_1_0_devc.csv', sep=',', header=1)
QP = M1.iloc[-1, 1]
ZP = M1.iloc[-1, 2]
errvec.append(np.abs(np.abs(QP) - vtexact))
Linf.append(np.abs(ZP - zexact(STIME)))

M2 = pd.read_csv(outdir + 'terminal_velocity_dt_0_1_devc.csv', sep=',', header=1)
QP = M2.iloc[-1, 1]
ZP = M2.iloc[-1, 2]
errvec.append(np.abs(np.abs(QP) - vtexact))
Linf.append(np.abs(ZP - zexact(STIME)))

M3 = pd.read_csv(outdir + 'terminal_velocity_dt_0_01_devc.csv', sep=',', header=1)
QP = M3.iloc[-1, 1]
ZP = M3.iloc[-1, 2]
errvec.append(np.abs(np.abs(QP) - vtexact))
Linf.append(np.abs(ZP - zexact(STIME)))

M4 = pd.read_csv(outdir + 'terminal_velocity_dt_0_001_devc.csv', sep=',', header=1)
QP = M4.iloc[-1, 1]
ZP = M4.iloc[-1, 2]
errvec.append(np.abs(np.abs(QP) - vtexact))
Linf.append(np.abs(ZP - zexact(STIME)))

M5 = pd.read_csv(outdir + 'terminal_velocity_dt_0_0001_devc.csv', sep=',', header=1)
QP = M5.iloc[-1, 1]
ZP = M5.iloc[-1, 2]
errvec.append(np.abs(np.abs(QP) - vtexact))
Linf.append(np.abs(ZP - zexact(STIME)))

# Convert to numpy arrays
errvec = np.array(errvec)
Linf = np.array(Linf)

if errvec[4] > 1e-6:
    print('Matlab Warning: The velocity in the terminal_velocity* cases is out of tolerance.')

if Linf[4] > 1e-3:
    print('Matlab Warning: The position in the terminal_velocity* cases is out of tolerance.')

fig = fdsplotlib.plot_to_fig(x_data=dtvec, y_data=errvec, marker_style='k*-', data_label='FDS',
                             x_min=1e-4, x_max=1, y_min=1e-15, y_max=1,
                             plot_type='loglog',
                             revision_label=version_string,
                             x_label='Time Step (s)',
                             y_label='Terminal Velocity Error (m/s)')

fdsplotlib.plot_to_fig(x_data=dtvec, y_data=dtvec,    marker_style='k--', data_label=r'$O(\delta t)$', figure_handle=fig)
fdsplotlib.plot_to_fig(x_data=dtvec, y_data=dtvec**2, marker_style='k-',  data_label=r'$O(\delta t^2)$', figure_handle=fig)

plt.savefig(pltdir + 'terminal_velocity_convergence.pdf', format='pdf')
plt.close()

# Create second figure - position convergence

fig = fdsplotlib.plot_to_fig(x_data=dtvec, y_data=Linf, marker_style='k*-', data_label='FDS',
                             x_min=1e-4, x_max=1, y_min=1e-10, y_max=1,
                             plot_type='loglog',
                             revision_label=version_string,
                             x_label='Time Step (s)',
                             y_label='Position Error (m/s)')

fdsplotlib.plot_to_fig(x_data=dtvec, y_data=dtvec,    marker_style='k--', data_label=r'$O(\delta t)$', figure_handle=fig)
fdsplotlib.plot_to_fig(x_data=dtvec, y_data=dtvec**2, marker_style='k-',  data_label=r'$O(\delta t^2)$', figure_handle=fig)

plt.savefig(pltdir + 'position_convergence.pdf', format='pdf')
plt.close()

