#!/usr/bin/env python
"""
McDermott
2-4-11
flame_height.py

integrates HRRPUL(z) from *_line.csv file to determine L_f/D (normalized flame height)
"""

import pandas as pd
import numpy as np
import os
import fdsplotlib

def check_hrr():
    # Get plot style parameters
    plot_style = fdsplotlib.get_plot_style('fds')

    outdir = '../../../out/Heskestad_Flame_Height/'
    box_file = os.path.join(outdir, 'box_height.csv')
    M = pd.read_csv(box_file, skiprows=1, header=None).values
    Qs = M[:, 0]
    Q = M[:, 1]
    RI = ['_RI=05', '_RI=10', '_RI=20']
    QI = ['p1','p2','p5','1','2','5','10','20','50','100','200','500','1000','2000','5000','10000']

    # Get Git version string if file exists
    git_file = os.path.join(outdir, 'Qs=1_RI=05_git.txt')
    version_string = fdsplotlib.get_version_string(git_file)

    # Create the figure
    fig = fdsplotlib.plot_to_fig(Qs, Q,
                                 x_min=0.05, x_max=1e4,
                                 y_min=1e2, y_max=1e8,
                                 x_label='$Q*$',
                                 y_label='Heat Release Rate (kW)',
                                 legend_location='lower right',
                                 plot_type='loglog',
                                 data_label='correct',
                                 plot_title='Flame Height Heat Release Verification',
                                 revision_label=version_string)

    # Add FDS results markers
    for j, q_val in enumerate(Qs):
        for i, ri in enumerate(RI):
            filename = f'Qs={QI[j]}{ri}_hrr.csv'
            filepath = os.path.join(outdir, filename)
            A = pd.read_csv(filepath, skiprows=2, header=None).values
            n = A.shape[0]
            Q_fds = np.mean(A[n//2:, 1])  # average over last half of rows

            marker_map = {0: 'ks', 1: 'r^', 2: 'go'}
            label_map = {0: r'$D^*/\delta x = 5$', 1: r'$D^*/\delta x = 10$', 2: r'$D^*/\delta x = 20$'}

            # Only add a data label for the first Qs entry
            if j == 0:
                label = label_map[i]
            else:
                label = None

            fdsplotlib.plot_to_fig([Qs[j]], [Q_fds],
                                   figure_handle=fig,
                                   marker_style=marker_map[i],
                                   data_label=label)

    # Save figure as PDF
    plotdir = '../../Manuals/FDS_Validation_Guide/SCRIPT_FIGURES/Heskestad/'
    fig_file = os.path.join(plotdir, 'Flame_Height_check_hrr.pdf')
    fig.savefig(fig_file, format='pdf')




# confirm heat release rate (uncomment or implement if needed)
check_hrr()

outdir = '../../../out/Heskestad_Flame_Height/'
expdir = '../../../exp/Heskestad_Flame_Height/'

filename = [
    ['Qs=p1_RI=05_devc.csv','Qs=p1_RI=10_devc.csv','Qs=p1_RI=20_devc.csv'],
    ['Qs=p2_RI=05_devc.csv','Qs=p2_RI=10_devc.csv','Qs=p2_RI=20_devc.csv'],
    ['Qs=p5_RI=05_devc.csv','Qs=p5_RI=10_devc.csv','Qs=p5_RI=20_devc.csv'],
    ['Qs=1_RI=05_devc.csv','Qs=1_RI=10_devc.csv','Qs=1_RI=20_devc.csv'],
    ['Qs=2_RI=05_devc.csv','Qs=2_RI=10_devc.csv','Qs=2_RI=20_devc.csv'],
    ['Qs=5_RI=05_devc.csv','Qs=5_RI=10_devc.csv','Qs=5_RI=20_devc.csv'],
    ['Qs=10_RI=05_devc.csv','Qs=10_RI=10_devc.csv','Qs=10_RI=20_devc.csv'],
    ['Qs=20_RI=05_devc.csv','Qs=20_RI=10_devc.csv','Qs=20_RI=20_devc.csv'],
    ['Qs=50_RI=05_devc.csv','Qs=50_RI=10_devc.csv','Qs=50_RI=20_devc.csv'],
    ['Qs=100_RI=05_devc.csv','Qs=100_RI=10_devc.csv','Qs=100_RI=20_devc.csv'],
    ['Qs=200_RI=05_devc.csv','Qs=200_RI=10_devc.csv','Qs=200_RI=20_devc.csv'],
    ['Qs=500_RI=05_devc.csv','Qs=500_RI=10_devc.csv','Qs=500_RI=20_devc.csv'],
    ['Qs=1000_RI=05_devc.csv','Qs=1000_RI=10_devc.csv','Qs=1000_RI=20_devc.csv'],
    ['Qs=2000_RI=05_devc.csv','Qs=2000_RI=10_devc.csv','Qs=2000_RI=20_devc.csv'],
    ['Qs=5000_RI=05_devc.csv','Qs=5000_RI=10_devc.csv','Qs=5000_RI=20_devc.csv'],
    ['Qs=10000_RI=05_devc.csv','Qs=10000_RI=10_devc.csv','Qs=10000_RI=20_devc.csv']
]

rho_inf = 1.2
cp = 1.0
T_inf = 293.0
g = 9.81
D = 1.13
f = 0.99

# Q* values from MATLAB comment
Qdot = [151, 303, 756, 1513, 3025, 7564, 15127, 30255,
        75636, 151273, 302545, 756363, 1512725, 3025450, 7563625, 15127250]

W = np.zeros((16,4))

# Loop through all heat release rates and resolutions
for i in range(16):  # hrr loop
    L = []
    for j in range(3):  # resolution loop
        df = pd.read_csv(os.path.join(outdir,filename[i][j]),skiprows=2,header=None)
        L.append(df.iloc[-1,3])  # 99th percentile (4th column)
        Qstar = Qdot[i]/(rho_inf*cp*T_inf*np.sqrt(g)*D**(5/2))
    W[i,:] = [Qstar,L[0]/D,L[1]/D,L[2]/D]

# Write file with FDS-predicted flame heights
header1 = ['Q*','L/D (RI=5)','L/D (RI=10)','L/D (RI=20)']
filename1 = os.path.join(outdir,'FDS_Flame_Height.csv')
pd.DataFrame(W,columns=header1).to_csv(filename1,index=False)

# Generate FDS results for Tamanini cases
filename_out = [
    os.path.join(outdir,'FDS_Tamanini_RI=05.csv'),
    os.path.join(outdir,'FDS_Tamanini_RI=10.csv'),
    os.path.join(outdir,'FDS_Tamanini_RI=20.csv')
]

fds_line_file = [
    ['Qs=1500_RI=05_line.csv','Qs=p6_RI=05_line.csv','Qs=p3_RI=05_line.csv'],
    ['Qs=1500_RI=10_line.csv','Qs=p6_RI=10_line.csv','Qs=p3_RI=10_line.csv'],
    ['Qs=1500_RI=20_line.csv','Qs=p6_RI=20_line.csv','Qs=p3_RI=20_line.csv']
]

Qstar_line = [1500,.6,.3]
header = ['z/L jet','Q jet','z/L 62','Q 62','z/L 31','Q 31']

for j in range(3):  # resolution loop
    A = None
    for k in range(3):  # hrr loop
        df = pd.read_csv(os.path.join(outdir,fds_line_file[j][k]),skiprows=2)
        z = df.iloc[:,0].values
        dz = z[1]-z[0]
        hrrpul = df.iloc[:,1].values
        Qdot_line = np.sum(hrrpul)*dz
        hrr = np.cumsum(hrrpul)*dz/Qdot_line
        kk = np.argmax(hrr>f)
        L = z[kk-1]+dz*(f-hrr[kk-1])/(hrr[kk]-hrr[kk-1])
        new_cols = np.column_stack(((z+dz/2)/L,hrr))
        A = new_cols if A is None else np.hstack((A,new_cols))

    pd.DataFrame(A,columns=header).to_csv(filename_out[j],index=False)

