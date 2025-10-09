
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.patches import FancyArrowPatch
import csv
import os
import fdsplotlib

plot_style = fdsplotlib.get_plot_style('fds')

outdir = '../../../out/Sandia_Plumes/'
expdir = '../../../exp/Sandia_Plumes/'
pltdir = '../../Manuals/FDS_Validation_Guide/SCRIPT_FIGURES/Sandia_Plumes/'

git_file = outdir + 'Sandia_CH4_1m_Test14_dx6cm_git.txt'
version_string = fdsplotlib.get_version_string(git_file)

devc_file = outdir + 'Sandia_CH4_1m_Test17_dx1p5cm_devc.csv'
tmin = 10
tmax = 20

M = np.genfromtxt(devc_file, delimiter=',', skip_header=2)

range_indices = np.where((M[:, 0] > tmin) & (M[:, 0] < tmax))[0]

t = M[range_indices, 0]
W = M[range_indices, 1]

fig = fdsplotlib.plot_to_fig(x_data=t, y_data=W, marker_style='k-',
                             x_min=10, x_max=20, y_min=0, y_max=7,
                             revision_label=version_string,
                             plot_title='Sandia Methane Pool Fire, Test 17',
                             x_label='Time (s)',
                             y_label='Vertical Velocity (m/s)')

fig.savefig(pltdir + 'Sandia_CH4_1m_Test17_dx1p5cm_velsignal.pdf', format='pdf')
plt.close(fig)


# Power spectrum

pave = None
devc_col = [1,2,3,4]

for i in devc_col:
    W = M[range_indices, i]
    n = len(W)
    u = np.fft.fft(W)
    p = u * np.conj(u) / n**2
    if i == devc_col[0]:
        pave = p.copy()
    else:
        pave = pave + p

pave = pave / len(devc_col)

dt_puff = 0.61
t = M[range_indices, 0]
T = t[len(t)-1] - t[0]
n2 = int(np.floor(n/2)) + 1
k = np.arange(2, n2+1) 
f = k / T
f_puff = 1 / dt_puff
kk = np.where((f > 5) & (f < 90))[0]

fig = fdsplotlib.plot_to_fig(x_data=f, y_data=pave[k-1].real, marker_style='k-',
                             x_min=0.1, x_max=1000, y_min=1e-6, y_max=1e1,
                             revision_label=version_string,
                             plot_title='Sandia Methane Pool Fire, Test 17',
                             plot_type='loglog',
                             x_label='Frequency (Hz)',
                             y_label='Autospectral Density (m$^2$/s)')

# Find peak frequency
k_fds = np.where(pave[k-1].real == np.max(pave[k-1].real))[0]
f_val = f[k_fds]

# Plot -5/3 slope
fdsplotlib.plot_to_fig(x_data=f[kk], y_data=f[kk]**(-5/3), marker_style='k-', figure_handle=fig)

# Plot puffing frequency line
fdsplotlib.plot_to_fig(x_data=[f_puff, f_puff], y_data=[np.min(pave[k-1].real), np.max(pave[k-1].real)], marker_style='k--', figure_handle=fig)

ax = plt.gca()
ax.text(0.15, 0.7, 'FDS W-Velocity, 1.5 cm Resolution', fontsize=plot_style['Label_Font_Size'])
ax.text(1.5e-1, 3.5e-4, 'measured', fontsize=plot_style['Label_Font_Size'])
ax.text(1.5e-1, 1.2e-4, 'puffing', fontsize=plot_style['Label_Font_Size'])
ax.text(1.5e-1, 0.4e-4, 'frequency', fontsize=plot_style['Label_Font_Size'])

arrow = FancyArrowPatch((0.3, 0.3), (0.4, 0.3), transform=fig.transFigure, arrowstyle='->', mutation_scale=20, linewidth=1.5, color='k')
fig.patches.append(arrow)
ax.text(1e1, 4e-2, '-5/3', fontsize=plot_style['Label_Font_Size'])

plt.savefig(pltdir + 'Sandia_CH4_1m_Test17_dx1p5cm_powerspectrum.pdf', format='pdf')
plt.close(fig)

