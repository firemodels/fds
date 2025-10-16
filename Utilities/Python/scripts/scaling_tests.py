
# Read _cpu.csv files for the MPI weak and strong scaling test cases

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import fdsplotlib

plot_style = fdsplotlib.get_plot_style('fds')

outdir = '../../../out/MPI_Scaling_Tests/'
pltdir = '../../Manuals/FDS_User_Guide/SCRIPT_FIGURES/'

git_file = outdir + 'strong_scaling_test_001_git.txt'
version_string = fdsplotlib.get_version_string(git_file)


M = []
M.append(pd.read_csv(outdir + 'strong_scaling_test_001_cpu.csv', skiprows=1, header=None))
M.append(pd.read_csv(outdir + 'strong_scaling_test_008_cpu.csv', skiprows=1, header=None))
M.append(pd.read_csv(outdir + 'strong_scaling_test_032_cpu.csv', skiprows=1, header=None))
M.append(pd.read_csv(outdir + 'strong_scaling_test_064_cpu.csv', skiprows=1, header=None))
M.append(pd.read_csv(outdir + 'strong_scaling_test_096_cpu.csv', skiprows=1, header=None))
M.append(pd.read_csv(outdir + 'strong_scaling_test_192_cpu.csv', skiprows=1, header=None))
M.append(pd.read_csv(outdir + 'strong_scaling_test_288_cpu.csv', skiprows=1, header=None))
M.append(pd.read_csv(outdir + 'strong_scaling_test_432_cpu.csv', skiprows=1, header=None))

r = np.array([1, 8, 32, 64, 96, 192, 288, 432])
r2 = np.array([0.1, 8, 32, 64, 96, 192, 432, 1000])

n_rows, n_cols = M[0].shape

# Initialize t and t2 arrays
t = np.zeros((8, n_cols))
t2 = np.zeros(8)

# Calculate normalized times
for j in range(n_cols):
    for i in range(8):
        t[i, j] = M[i].iloc[0, j] / M[0].iloc[0, n_cols - 1]
        t2[i] = 1.0 / r2[i]

# Check if strong scaling test is within tolerance
if t[7, n_cols - 1] > 4 / r[7] or t[7, n_cols - 1] < 1 / r[7]:
    print('Error: strong scaling test out of tolerance')

fig = fdsplotlib.plot_to_fig(x_data=[1e-10,1e-10], y_data=[1e-10,1e-10],
                             x_min=1, x_max=1000, y_min=0.00001, y_max=1,
                             plot_type='loglog',
                             legend_location='outside',
                             revision_label=version_string,
                             plot_title='Strong Scaling Test',
                             x_label='MPI Processes',
                             y_label='Relative Wall Clock Time')

for i in range(10):
    fdsplotlib.plot_to_fig(x_data=r2, y_data=t2/2**(i-1), figure_handle=fig, marker_style='k:') 

marks = np.array(['k-o','r-o','b-o','m-o','c-o','g-o','y-o','k-s'])
ind = np.array([n_cols-1,2,3,4,5,11,9,1])
legend = np.array(['Total', 'DIVG', 'MASS', 'VELO', 'PRES', 'COMM', 'RADI', 'MAIN'])

for i in range(8):
    fdsplotlib.plot_to_fig(x_data=r, y_data=t[:,ind[i]], figure_handle=fig, marker_style=marks[i], data_label=legend[i]) 

plt.savefig(pltdir + 'strong_scaling_test.pdf', format='pdf')
plt.close()

# Import data files for weak scaling test

M = []
M.append(pd.read_csv(outdir + 'weak_scaling_test_001_cpu.csv', skiprows=1, header=None))
M.append(pd.read_csv(outdir + 'weak_scaling_test_002_cpu.csv', skiprows=1, header=None))
M.append(pd.read_csv(outdir + 'weak_scaling_test_004_cpu.csv', skiprows=1, header=None))
M.append(pd.read_csv(outdir + 'weak_scaling_test_008_cpu.csv', skiprows=1, header=None))
M.append(pd.read_csv(outdir + 'weak_scaling_test_016_cpu.csv', skiprows=1, header=None))
M.append(pd.read_csv(outdir + 'weak_scaling_test_032_cpu.csv', skiprows=1, header=None))
M.append(pd.read_csv(outdir + 'weak_scaling_test_064_cpu.csv', skiprows=1, header=None))
M.append(pd.read_csv(outdir + 'weak_scaling_test_128_cpu.csv', skiprows=1, header=None))
M.append(pd.read_csv(outdir + 'weak_scaling_test_192_cpu.csv', skiprows=1, header=None))
M.append(pd.read_csv(outdir + 'weak_scaling_test_288_cpu.csv', skiprows=1, header=None))
M.append(pd.read_csv(outdir + 'weak_scaling_test_432_cpu.csv', skiprows=1, header=None))

r = np.array([1, 2, 4, 8, 16, 32, 64, 128, 192, 288, 432])

n_rows, n_cols = M[0].shape

t = np.zeros(11)
t2 = np.ones(11)

for i in range(11):
    t[i] = M[0].iloc[0, n_cols - 1] / M[i].iloc[0, n_cols - 1]

if t[10] > 1.0 or t[10] < 0.4:
    print('Error: weak scaling test out of tolerance')

fig = fdsplotlib.plot_to_fig(x_data=r, y_data=t, data_label='FDS', marker_style='ko',
                             x_min=1, x_max=1000, y_min=0, y_max=1.2,
                             plot_type='semilogx',
                             revision_label=version_string,
                             plot_title='Weak Scaling Test',
                             x_label='MPI Processes',
                             y_label='Efficiency')

fdsplotlib.plot_to_fig(x_data=r, y_data=t2, data_label='Ideal', marker_style='k--', figure_handle=fig)

plt.savefig(pltdir + 'weak_scaling_test.pdf', format='pdf')
plt.close()

