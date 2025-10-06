
# Contour plots of vertical centerline temperature and back-layer location vs time for Memorial Tunnel cases with transverse (i.e. no fans) ventilation.

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.interpolate import RegularGridInterpolator
import os
import sys
import fdsplotlib

plot_style = fdsplotlib.get_plot_style('fds')

outdir = '../../../out/Memorial_Tunnel/'
expdir = '../../../exp/Memorial_Tunnel/'
pltdir = '../../Manuals/FDS_Validation_Guide/SCRIPT_FIGURES/Memorial_Tunnel/'

git_file = outdir + 'Test_102_cat_git.txt'
version_string = fdsplotlib.get_version_string(git_file)

times = [10]  # Times (min) to make the centerline temperature plots

# Loop  214   213    211    209    208    207    307    306    305    205    304    303    302    301    202
pos  = [19.8, 125.6, 230.7, 340.8, 446.2, 528.2, 573.3, 586.1, 604.1, 615.4, 627.6, 645.0, 681.5, 723.3, 833.9]
hgt_mod = [[0.3, 1.1, 1.8, 2.4, 3.0, 3.5, 4.0], [0.3, 1.2, 2.4, 3.7, 4.8, 5.7, 6.5, 7.0, 7.4]]
hgt_exp = [[0.3, 1.1, 1.8, 2.4, 3.0, 3.5, 4.0], [0.3, 1.2, 2.4, 3.7, 4.8, 5.7, 6.5, 7.0]]

test  = ['101CR', '102', '102R1', '102R', '103', '104', '105', '106', '107', '108',
         '109', '110', '111', '112A', '113A', '115A', '126B', '126BR1', '128B', '202',
         '203', '205', '207A', '208A', '210', '212', '214A', '215A', '216A', '217A',
         '218B', '223', '226', '227A', '229', '230', '231', '233', '235', '236',
         '238A', '239', '244B', '245B', '246B', '247B', '248B', '249B', '250B', '251B',
         '252B', '301A', '302A', '303A', '305A', '306A', '309A', '312A', '313A', '314',
         '315A', '316', '317A', '318A', '319A', '320A', '321A', '338B', '339B', '340B',
         '341B', '342B', '343B', '344B', '345B', '346B', '401A', '403A', '404A', '407B',
         '408B']

sequence = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
            1, 1, 1, 1, 1, 1, 1, 1, 1, 3,
            3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
            3, 4, 5, 5, 6, 6, 6, 6, 6, 6,
            6, 6, 5, 6, 6, 3, 3, 3, 6, 6,
            3, 8, 8, 8, 8, 8, 8, 8, 8, 9,
            101, 9, 102, 101, 102, 102, 9, 8, 8, 8,
            8, 8, 8, 8, 8, 8, 13, 13, 14, 13,
            13]

levels = [50, 100, 200, 400, 600, 800]  # Temperatures contours (C)
single_level = [50]                     # 50 C is the indicator of a back-layer

# T Loops: 214, 213, 211, 209, 208, 207, 307, 306, 305, 205, 304, 303, 302, 301, 202
mod_data_indices = [list(range(1,8)), list(range(15,22)), list(range(22,29)), list(range(29,36)), list(range(43,50)), list(range(57,64)), list(range(71,78)), list(range(85,92)), list(range(92,99)), list(range(106,113)), list(range(113,120)), list(range(127,134)), list(range(134,141)), list(range(148,155)), list(range(162,169))]
exp_data_indices = [list(range(99,106)), list(range(92,99)), list(range(85,92)), list(range(78,85)), list(range(71,78)), list(range(64,71)), list(range(57,64)), list(range(50,57)), list(range(43,50)), list(range(36,43)), list(range(29,36)), list(range(22,29)), list(range(15,22)), list(range(8,15)), list(range(1,8))]

def read_csv_with_header(filename):
    """
    Read CSV file similar to Matlab's importdata
    Returns object with data attribute and colheaders attribute
    """
    class ImportData:
        def __init__(self):
            self.data = None
            self.colheaders = None

    result = ImportData()
    try:
        # Read first two lines as header
        df = pd.read_csv(filename, skiprows=2, header=None)
        result.data = df.values
        # Read column headers
        with open(filename, 'r') as f:
            lines = f.readlines()
            if len(lines) >= 2:
                result.colheaders = lines[1].strip().split(',')
    except:
        return None

    return result

# Make contour plots of centerline tunnel temperature for all experiments and 19 individual times.
# Also, make a contour plot for each experiment showing a single temperature contour as a function of position and time.

for k in range(81):  # Experiments

    M = None
    E = None
    H = None

    try:
        M  = read_csv_with_header(os.path.join(outdir, f'Test_{test[k]}_cat_devc.csv'))
        E  = read_csv_with_header(os.path.join(expdir, f'TP{test[k]}.csv'))
        H  = read_csv_with_header(os.path.join(expdir, f'HR{test[k]}.csv'))
    except:
        continue

    if M is None or E is None or H is None:
        continue

    if M.data[-1, 0] < 200:
        continue

    if sequence[k] == 1:
        ventilation_type = 'Full Transverse Ventilation'
    elif sequence[k] == 3:
        ventilation_type = 'Partial Transverse Exhaust Ventilation'
    elif sequence[k] in [4, 5]:
        ventilation_type = 'Partial Transverse Exhaust Ventilation'
    elif sequence[k] == 6:
        ventilation_type = 'Two-Zone Partial Transverse Ventilation'
    elif sequence[k] == 8:
        ventilation_type = 'Partial Transverse Ventilation with Single Point Extraction'
    elif sequence[k] == 9:
        ventilation_type = 'Point Supply Operation'
    elif sequence[k] == 10:
        ventilation_type = 'Point Exhaust Operation'
    elif sequence[k] in [13, 14]:
        ventilation_type = 'Partial Transverse Ventilation with Oversize Exhaust Ports'
    else:
        ventilation_type = ''

    for i in range(len(times)):  # Times

        X_mod = None
        Y_mod = None
        Z_mod = None
        X_exp = None
        Y_exp = None
        Z_exp = None

        # Get time indices using interpolation
        mod_time_index = int(np.interp(60*times[i], M.data[:, 0], np.arange(len(M.data[:, 0]))))
        exp_time_index = int(np.interp(60*times[i], E.data[:, 0], np.arange(len(E.data[:, 0]))))
        hrr_time_index = int(np.interp(60*times[i], H.data[:, 0], np.arange(len(H.data[:, 0]))))

        X_mod, Y_mod = np.meshgrid(pos[0:15], hgt_mod[0])
        X_exp, Y_exp = np.meshgrid(pos[0:15], hgt_exp[0])

        Z_mod = np.zeros((len(hgt_mod[0]), 15))
        Z_exp = np.zeros((len(hgt_exp[0]), 15))

        for kk in range(15):  # Loops
            mod_indices = [idx for idx in mod_data_indices[kk]]
            exp_indices = [idx for idx in exp_data_indices[kk]]
            exp_indices_flipped = list(reversed(exp_indices))

            Z_mod[:, kk] = M.data[mod_time_index, mod_indices]
            Z_exp[:, kk] = E.data[exp_time_index, exp_indices_flipped]

        newpoints = 100
        X_mod_interp, Y_mod_interp = np.meshgrid(
            np.linspace(X_mod.min(), X_mod.max(), newpoints),
            np.linspace(Y_mod.min(), Y_mod.max(), newpoints)
        )

        f_mod = RegularGridInterpolator( (Y_mod[:, 0], X_mod[0, :]), Z_mod)
        Z_mod_interp = f_mod((Y_mod_interp, X_mod_interp))

        X_exp_interp, Y_exp_interp = np.meshgrid(
            np.linspace(X_exp.min(), X_exp.max(), newpoints),
            np.linspace(Y_exp.min(), Y_exp.max(), newpoints)
        )

        f_exp = RegularGridInterpolator( (Y_exp[:, 0], X_exp[0, :]), Z_exp) 
        Z_exp_interp = f_exp((Y_exp_interp, X_exp_interp))

        fig = fdsplotlib.plot_to_fig(x_data=[-1,-1], y_data=[-1,-1],
                                     x_min=0, x_max=854, y_min=0, y_max=4.4,
                                     figure_size=(6.5,1.5),
                                     plot_size=(6.0,1.0),
                                     plot_origin=(0.4,0.35),
                                     axeslabel_fontsize=7,
                                     ticklabel_fontsize=7,
                                     title_fontsize=10,
                                     revision_label=version_string,
                                     version_fontsize=7,
                                     x_label='Tunnel Length (m)',
                                     y_label='Height (m)')

        CS_mod = plt.contour(X_mod_interp, Y_mod_interp, Z_mod_interp, levels=levels, colors='red', linewidths=1)
        CS_exp = plt.contour(X_exp_interp, Y_exp_interp, Z_exp_interp, levels=levels, colors='black', linewidths=1)
        plt.clabel(CS_mod, fontsize=6, colors='red')
        plt.clabel(CS_exp, fontsize=6, colors='black')

        ax = plt.gca()
        ax.set_xticks([0, 100, 200, 300, 400, 500, 600, 700, 800, 854])
        ax.set_xticklabels(['0', '100', '200', '300', '400', '500', '600', '700', '800', '854'])
        ax.set_yticks([0,1,2,3,4])

        ax.text(10, 4.0, f'Test {test[k]}', fontsize=7)
        ax.text(10, 3.5, f'Time: {times[i]} min', fontsize=7)
        ax.text(10, 3.0, f'HRR: {H.data[hrr_time_index, 1]/1000.:.1f} MW', fontsize=7)

        from matplotlib.lines import Line2D
        legend_elements = [ Line2D([0], [0], color='red', label='FDS', linewidth=1),
                            Line2D([0], [0], color='black', label='Exp', linewidth=1) ]
        plt.legend(handles=legend_elements, loc='lower left', fontsize=7)

        plt.savefig(os.path.join(pltdir, f'Test_{test[k]}_T_{times[i]}.pdf'), format='pdf')
        plt.close(fig)

    # For each experiment, make a contour plot of the extent of a single temperature contour at each time during the experiment

    X_mod = None
    Y_mod = None
    Z_mod = None
    X_exp = None
    Y_exp = None
    Z_exp = None

    X_mod, Y_mod = np.meshgrid(pos[0:15], M.data[:, 0]/60)
    X_exp, Y_exp = np.meshgrid(pos[0:15], E.data[:, 0]/60)

    Z_mod = np.zeros((len(M.data[:, 0]), 15))
    Z_exp = np.zeros((len(E.data[:, 0]), 15))

    for kk in range(len(M.data[:, 0])):
        for ii in range(15):
            mod_idx = mod_data_indices[ii][6] 
            Z_mod[kk, ii] = M.data[kk, mod_idx]

    for kk in range(len(E.data[:, 0])):
        for ii in range(15):
            exp_idx = exp_data_indices[ii][0]
            Z_exp[kk, ii] = E.data[kk, exp_idx]

    newpoints = 100
    X_mod_interp, Y_mod_interp = np.meshgrid(
        np.linspace(X_mod.min(), X_mod.max(), newpoints),
        np.linspace(Y_mod.min(), Y_mod.max(), newpoints)
    )

    f_mod = RegularGridInterpolator( (Y_mod[:, 0], X_mod[0, :]), Z_mod) 
    Z_mod_interp = f_mod((Y_mod_interp, X_mod_interp))

    X_exp_interp, Y_exp_interp = np.meshgrid(
        np.linspace(X_exp.min(), X_exp.max(), newpoints),
        np.linspace(Y_exp.min(), Y_exp.max(), newpoints)
    )

    f_exp = RegularGridInterpolator( (Y_exp[:, 0], X_exp[0, :]), Z_exp) 
    Z_exp_interp = f_exp((Y_exp_interp, X_exp_interp))

    fig = fdsplotlib.plot_to_fig(x_data=[615.4,615.4], y_data=[0,30], marker_style='k--',
                                 x_min=0, x_max=854, y_min=0, y_max=30,
                                 figure_size=(6.5,1.5),
                                 plot_size=(6.0,1.0),
                                 plot_origin=(0.4,0.35),
                                 axeslabel_fontsize=7,
                                 ticklabel_fontsize=7,
                                 title_fontsize=10,
                                 revision_label=version_string,
                                 version_fontsize=7,
                                 x_label='Tunnel Length (m)',
                                 y_label='Time (min)')

    CS_mod = plt.contour(X_mod_interp, Y_mod_interp, Z_mod_interp, levels=single_level, colors='red', linewidths=1)
    CS_exp = plt.contour(X_exp_interp, Y_exp_interp, Z_exp_interp, levels=single_level, colors='black', linewidths=1)
    plt.clabel(CS_mod, fontsize=6, colors='red')
    plt.clabel(CS_exp, fontsize=6, colors='black')

    ax = plt.gca()
    ax.set_xticks([0, 100, 200, 300, 400, 500, 600, 700, 800, 854])
    ax.set_xticklabels(['0', '100', '200', '300', '400', '500', '600', '700', '800', '854'])
    ax.set_yticks([0, 10, 20, 30])

    ax.text(10, 26, f'Test {test[k]}', fontsize=7)
    ax.text(10, 22, ventilation_type, fontsize=7)

    from matplotlib.lines import Line2D
    legend_elements = [ Line2D([0], [0], color='red', label='FDS', linewidth=1),
                        Line2D([0], [0], color='black', label='Exp', linewidth=1) ]
    plt.legend(handles=legend_elements, loc='lower left', fontsize=7)

    plt.savefig(os.path.join(pltdir, f'Test_{test[k]}_tvT.pdf'), format='pdf')
    plt.close(fig)

