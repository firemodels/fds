
# Show efficiency of OpenMP timing cases

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import fdsplotlib

plot_style = fdsplotlib.get_plot_style('fds')

outdir = '../../../out/OMP_Scaling_Tests/'
pltdir = '../../Manuals/FDS_User_Guide/SCRIPT_FIGURES/'

git_file = outdir + 'openmp_test64a_git.txt'
version_string = fdsplotlib.get_version_string(git_file)

ncores = [1, 2, 3, 4, 5, 6, 7, 8]
a = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h']

time64 = np.zeros(len(a))
time128 = np.zeros(len(a))

# Read data for 64 cube cases
for i in range(len(a)):
    M = pd.read_csv(outdir + f'openmp_test64{a[i]}_devc.csv', skiprows=1)
    M.columns = M.columns.str.strip()
    j = M.columns.get_loc('clock time')
    time64[i] = M.iloc[-1, j]
time64 = time64 / time64[0] * 100

for i in range(len(a)):
    M = pd.read_csv(outdir + f'openmp_test128{a[i]}_devc.csv', skiprows=1)
    M.columns = M.columns.str.strip()
    j = M.columns.get_loc('clock time')
    time128[i] = M.iloc[-1, j]
time128 = time128 / time128[0] * 100

fig = fdsplotlib.plot_to_fig(x_data=ncores, y_data=time64, marker_style='b^-', data_label='$64^3$',
                             x_min=1, x_max=8, y_min=0, y_max=100,
                             revision_label=version_string,
                             x_label='Number of OpenMP Threads',
                             y_label=r'Relative Clock Time (%)')

fdsplotlib.plot_to_fig(x_data=ncores, y_data=time128, marker_style='rsq-', data_label='$128^3$', figure_handle=fig)

plt.savefig(pltdir + 'openmp_timing_benchmarks.pdf', format='pdf')
plt.close()

if time64[3] > 80.:
    print(f'Matlab Warning: Timing for openmp_test64 out of tolerance. {time64[3]}')

if time128[3] > 80.:
    print(f'Matlab Warning: Timing for openmp_test128 out of tolerance. {time128[3]}')

