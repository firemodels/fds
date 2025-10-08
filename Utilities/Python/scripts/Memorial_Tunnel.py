
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from scipy.interpolate import interp1d, griddata
import os
import fdsplotlib

plot_style = fdsplotlib.get_plot_style('fds')

expdir = '../../../exp/Memorial_Tunnel/'
outdir = '../../../out/Memorial_Tunnel/'
pltdir = '../../Manuals/FDS_Validation_Guide/SCRIPT_FIGURES/Memorial_Tunnel/'

git_file = outdir + 'Cold_Flow_Series_1_cat_git.txt'
version_string = fdsplotlib.get_version_string(git_file)


def importdata(filename, delimiter=',', skip_header=2):
    """
    Import data from CSV file similar to Matlab's importdata

    Parameters:
    -----------
    filename : str
        Path to the file
    delimiter : str
        Delimiter character
    skip_header : int
        Number of header rows to skip

    Returns:
    --------
    Object with 'data' attribute containing numpy array
    """
    class ImportData:
        def __init__(self):
            self.data = None

    result = ImportData()

    try:
        result.data = np.genfromtxt(filename, delimiter=delimiter, skip_header=skip_header)
    except Exception as e:
        print(f"Error loading {filename}: {e}")
        result.data = np.array([[]])

    return result

times = [[5, 10, 16, 20, 26, 30],      # 501
         [5, 10, 16, 20, 26, 30],      # 502
         [12, 16, 22, 26, 28, 30],     # 605
         [4, 8, 16, 20, 24, 26],       # 606A
         [8, 10, 14, 22, 26, 28],      # 607
         [1, 2, 4, 6, 12, 20],         # 608
         [6, 8, 12, 22, 24, 26],       # 610
         [2, 6, 8, 16, 20, 26],        # 611
         [2, 6, 10, 14, 18, 22],       # 612B
         [2, 4, 6, 10, 18, 28],        # 615B
         [16, 18, 20, 22, 28, 30],     # 617A
         [2, 5, 6, 12, 24, 26],        # 618A
         [1, 2, 3, 12, 14, 24],        # 621A
         [3, 4, 5, 8, 10, 14],         # 622B
         [2, 4, 6, 10, 20, 28],        # 623B
         [1, 2, 3, 4, 6, 8],           # 624B
         [6, 8, 10, 18, 20, 21]]       # 625B

pos = [19.8, 125.6, 230.7, 340.8, 446.2, 528.2, 573.3, 586.1, 604.1, 615.4, 627.6, 645.0, 681.5, 723.3, 833.9]
hgt_mod = [[0.3, 1.1, 2.0, 2.6, 3.2, 3.7, 4.1], [0.3, 1.2, 2.4, 3.7, 4.8, 5.7, 6.5, 7.0, 7.4]]
hgt_exp = [[0.3, 1.1, 2.0, 2.6, 3.2, 3.7, 4.1], [0.3, 1.2, 2.4, 3.7, 4.8, 5.7, 6.5, 7.0]]
test = ['501', '502', '605', '606A', '607', '608', '610', '611', '612B', '615B', '617A', '618A', '621A', '622B', '623B', '624B', '625B']
levels = [50, 100, 200, 400, 600, 800]
single_level = [50]

# Loops:             214         , 213         , 211        , 209       , 208       , 207        , 307         , 306         , 305         , 205         , 304         , 303         , 302         , 301         , 202
mod_data_indices = [list(range(1, 8)), list(range(15, 24)), list(range(24, 33)), list(range(33, 42)), list(range(51, 60)), list(range(69, 78)), list(range(87, 96)), list(range(105, 114)), list(range(114, 123)), list(range(132, 141)), list(range(141, 150)), list(range(159, 168)), list(range(168, 177)), list(range(186, 195)), list(range(204, 213))]
exp_data_indices = [list(range(112, 119)), list(range(104, 112)), list(range(96, 104)), list(range(88, 96)), list(range(80, 88)), list(range(72, 80)), list(range(64, 72)), list(range(56, 64)), list(range(48, 56)), list(range(40, 48)), list(range(32, 40)), list(range(24, 32)), list(range(16, 24)), list(range(8, 16)), list(range(1, 8))]

# Make contour plots of centerline tunnel temperature for 17 experiments and 19 individual times.
# Also, make a contour plot for each experiment showing a single temperature contour as a function of position and time.

for k in range(17):  # Experiments

    M = importdata(outdir + 'Test_' + test[k] + '_cat_devc.csv', ',', 2)
    E = importdata(expdir + 'TP' + test[k] + '.csv', ',', 2)
    EV = importdata(expdir + 'QP' + test[k] + '.csv', ',', 2)
    H = importdata(expdir + 'HRR' + test[k] + '.csv', ',', 2)

    for i in range(len(times[k])):  # Times

        # Interpolate time indices
        mod_time_index = int(np.interp(60 * times[k][i], M.data[:, 0], np.arange(len(M.data[:, 0]))))
        exp_time_index = int(np.interp(60 * times[k][i], E.data[:, 0], np.arange(len(E.data[:, 0]))))
        hrr_time_index = int(np.interp(60 * times[k][i], H.data[:, 0], np.arange(len(H.data[:, 0]))))
        exp_VF_time_index = int(np.interp(60 * times[k][i], EV.data[:, 0], np.arange(len(EV.data[:, 0]))))

        fig = fdsplotlib.plot_to_fig(x_data=[-1,-1], y_data=[-1,-1],
                                     x_min=0, x_max=854, y_min=0, y_max=8,
                                     figure_size=(6.5,1.5),
                                     plot_size=(6.0,1.0),
                                     plot_origin=(0.4,0.3),
                                     axeslabel_fontsize=7,
                                     ticklabel_fontsize=7,
                                     title_fontsize=10,
                                     revision_label=version_string,
                                     version_fontsize=7,
                                     x_label='Tunnel Length (m)',
                                     y_label='Height (m)')

        # Create meshgrids for model and experiment data
        X_mod, Y_mod = np.meshgrid(pos[1:14], hgt_mod[1])
        X_exp, Y_exp = np.meshgrid(pos[1:14], hgt_exp[1])

        # Extract data for each loop
        Z_mod = np.zeros((len(hgt_mod[1]), 13))
        Z_exp = np.zeros((len(hgt_exp[1]), 13))

        for kk in range(1, 14):  # Loops (Matlab: 2:14, Python: 1:14)
            Z_mod[:, kk - 1] = M.data[mod_time_index, mod_data_indices[kk]]
            Z_exp[:, kk - 1] = E.data[exp_time_index, exp_data_indices[kk][::-1]]

        # Interpolate data
        newpoints = 100
        X_mod_interp, Y_mod_interp = np.meshgrid(
            np.linspace(np.min(X_mod), np.max(X_mod), newpoints),
            np.linspace(np.min(Y_mod), np.max(Y_mod), newpoints)
        )

        # Flatten arrays for griddata interpolation
        points_mod = np.array([X_mod.flatten(), Y_mod.flatten()]).T
        values_mod = Z_mod.flatten()
        Z_mod_interp = griddata(points_mod, values_mod, (X_mod_interp, Y_mod_interp), method='cubic')

        X_exp_interp, Y_exp_interp = np.meshgrid(
            np.linspace(np.min(X_exp), np.max(X_exp), newpoints),
            np.linspace(np.min(Y_exp), np.max(Y_exp), newpoints)
        )

        points_exp = np.array([X_exp.flatten(), Y_exp.flatten()]).T
        values_exp = Z_exp.flatten()
        Z_exp_interp = griddata(points_exp, values_exp, (X_exp_interp, Y_exp_interp), method='cubic')

        C_mod = plt.contour(X_mod_interp, Y_mod_interp, Z_mod_interp, levels, colors='red', linestyles='solid', linewidths=1)
        C_exp = plt.contour(X_exp_interp, Y_exp_interp, Z_exp_interp, levels, colors='black', linestyles='solid', linewidths=1)
        plt.clabel(C_mod, fontsize=6, colors='red')
        plt.clabel(C_exp, fontsize=6, colors='black')

        ax = plt.gca()
        ax.tick_params(axis='both', which='major', labelsize=7, length=0)
        ax.set_xticks([0, 100, 200, 300, 400, 500, 600, 700, 800, 854])
        ax.set_xticklabels(['0', '100', '200', '300', '400', '500', '600', '700', '800', '854'])

        ax.text(10, 7.2, 'Test ' + test[k], fontsize=7)
        ax.text(10, 6.3, 'Time: ' + str(times[k][i]) + ' min', fontsize=7)
        ax.text(10, 5.4, 'HRR: ' + '{:.1f}'.format(H.data[hrr_time_index, 1] / 1000.) + ' MW', fontsize=7)

        from matplotlib.lines import Line2D
        legend_elements = [ Line2D([0], [0], color='red', label='FDS', linewidth=1),
                            Line2D([0], [0], color='black', label='Exp', linewidth=1) ]
        plt.legend(handles=legend_elements, loc='upper right', fontsize=7)

        plt.savefig(pltdir + 'Test_' + test[k] + '_T_' + str(times[k][i]) + '.pdf', format='pdf')
        plt.close()

    # For each experiment, make a contour plot of the extent of a single temperature contour at each time during the experiment

    X_mod, Y_mod = np.meshgrid(pos[1:14], M.data[:, 0] / 60)
    X_exp, Y_exp = np.meshgrid(pos[1:14], E.data[:, 0] / 60)

    Z_mod = np.zeros((len(M.data[:, 0]), 13))
    Z_exp = np.zeros((len(E.data[:, 0]), 13))

    for kk in range(len(M.data[:, 0])):
        for ii in range(1, 14):
            Z_mod[kk, ii - 1] = M.data[kk, mod_data_indices[ii][7]]

    for kk in range(len(E.data[:, 0])):
        for ii in range(1, 14):
            Z_exp[kk, ii - 1] = E.data[kk, exp_data_indices[ii][0]]

    # Interpolate data
    newpoints = 100
    X_mod_interp, Y_mod_interp = np.meshgrid(
        np.linspace(np.min(X_mod), np.max(X_mod), newpoints),
        np.linspace(np.min(Y_mod), np.max(Y_mod), newpoints)
    )

    points_mod = np.array([X_mod.flatten(), Y_mod.flatten()]).T
    values_mod = Z_mod.flatten()
    Z_mod_interp = griddata(points_mod, values_mod, (X_mod_interp, Y_mod_interp), method='cubic')

    X_exp_interp, Y_exp_interp = np.meshgrid(
        np.linspace(np.min(X_exp), np.max(X_exp), newpoints),
        np.linspace(np.min(Y_exp), np.max(Y_exp), newpoints)
    )

    points_exp = np.array([X_exp.flatten(), Y_exp.flatten()]).T
    values_exp = Z_exp.flatten()
    Z_exp_interp = griddata(points_exp, values_exp, (X_exp_interp, Y_exp_interp), method='cubic')

    fig = fdsplotlib.plot_to_fig(x_data=[615.4,615.4], y_data=[0,30], marker_style='k--',
                                 x_min=0, x_max=854, y_min=0, y_max=30,
                                 figure_size=(6.5,1.3),
                                 plot_size=(6.0,1.0),
                                 plot_origin=(0.4,0.1),
                                 axeslabel_fontsize=7,
                                 ticklabel_fontsize=7,
                                 title_fontsize=10,
                                 revision_label=version_string,
                                 version_fontsize=7,
                                 y_label='Time (min)')

    C_mod = plt.contour(X_mod_interp, Y_mod_interp, Z_mod_interp, single_level, colors='red', linestyles='solid', linewidths=1)
    C_exp = plt.contour(X_exp_interp, Y_exp_interp, Z_exp_interp, single_level, colors='black', linestyles='solid', linewidths=1)
    plt.clabel(C_mod, fontsize=6, colors='red')
    plt.clabel(C_exp, fontsize=6, colors='black')

    ax = plt.gca()
    ax.tick_params(axis='both', which='major', labelsize=7, length=0)
    ax.set_xticks(pos)
    ax.set_xticklabels(['', '', '', '', '', '', '', '', '', '', '', '', '', '', ''])

    ax.text(10, 26, 'Test ' + test[k], fontsize=7)
    #ax.text(10, 22, 'FDS red; Exp black', fontsize=7)

    from matplotlib.lines import Line2D
    legend_elements = [ Line2D([0], [0], color='red', label='FDS', linewidth=1),
                        Line2D([0], [0], color='black', label='Exp', linewidth=1) ]
    plt.legend(handles=legend_elements, loc='upper right', fontsize=7)

    plt.savefig(pltdir + 'Test_' + test[k] + '_tvT.pdf', format='pdf')
    plt.close()

# End Experiment loop


# Collect critical velocity information

cv_time = [[10, 12, 14, 16, 18, 20],  # 501
           [10, 12, 14, 16, 18, 20],  # 502
           [25, 28],                  # 605
           [24],                      # 606A
           [2, 10],                   # 607
           [],                        # 608
           [2, 8, 20, 25],            # 610
           [10, 16],                  # 611
           [15, 20],                  # 612B
           [13, 20],                  # 615B
           [14, 17, 27, 30],          # 617A
           [20, 26],                  # 618A
           [8, 12],                   # 621A
           [3, 6],                    # 622B
           [],                        # 623B
           [10, 14],                  # 624B
           []]                        # 625B

js = 0
jn = 0
ks = 0
kn = 0

V_crit_exp_sup = []
HRR_crit_exp_sup = []
V_crit_exp_nosup = []
HRR_crit_exp_nosup = []
V_crit_mod_sup = []
HRR_crit_mod_sup = []
V_crit_mod_nosup = []
HRR_crit_mod_nosup = []

for k in range(2, 17):  # Experiments (Matlab: 3:17, Python: 2:17)

    M = importdata(outdir + 'Test_' + test[k] + '_cat_devc.csv', ',', 2)
    E = importdata(expdir + 'TP' + test[k] + '.csv', ',', 2)
    EV = importdata(expdir + 'QP' + test[k] + '.csv', ',', 2)
    H = importdata(expdir + 'HRR' + test[k] + '.csv', ',', 2)
    HM = importdata(outdir + 'Test_' + test[k] + '_cat_hrr.csv', ',', 2)

    for i in range(len(cv_time[k])):  # Times

        mod_time_index = int(np.interp(60 * cv_time[k][i], M.data[:, 0], np.arange(len(M.data[:, 0]))))
        exp_time_index = int(np.interp(60 * cv_time[k][i], E.data[:, 0], np.arange(len(E.data[:, 0]))))
        hrr_time_index = int(np.interp(60 * cv_time[k][i], H.data[:, 0], np.arange(len(H.data[:, 0]))))
        hrr_mod_time_index = int(np.interp(60 * cv_time[k][i], HM.data[:, 0], np.arange(len(HM.data[:, 0]))))
        exp_VF_time_index = int(np.interp(60 * cv_time[k][i], EV.data[:, 0], np.arange(len(EV.data[:, 0]))))

        if E.data[exp_time_index, 63] < 30:
            V_crit_exp_sup.append(-EV.data[exp_VF_time_index, 1] / 60.4)
            HRR_crit_exp_sup.append(H.data[hrr_time_index, 1] / 1000)
            js = js + 1
        else:
            V_crit_exp_nosup.append(-EV.data[exp_VF_time_index, 1] / 60.4)
            HRR_crit_exp_nosup.append(H.data[hrr_time_index, 1] / 1000)
            jn = jn + 1

        if M.data[mod_time_index, 112] < 30:
            V_crit_mod_sup.append(-M.data[mod_time_index, 218] / 60.4)
            HRR_crit_mod_sup.append(HM.data[hrr_mod_time_index, 1] / 1000)
            ks = ks + 1
        else:
            V_crit_mod_nosup.append(-M.data[mod_time_index, 218] / 60.4)
            HRR_crit_mod_nosup.append(HM.data[hrr_mod_time_index, 1] / 1000)
            kn = kn + 1

# Create critical velocity plot

fig = fdsplotlib.plot_to_fig(x_data=[8.5,105], y_data=[1.9,3.5], marker_style='k-',
                             data_label='Theoretical Critical Velocity',
                             x_min=5, x_max=200, y_min=0, y_max=4,
                             plot_type='semilogx',
                             revision_label=version_string,
                             x_label='Heat Release Rate (MW)',
                             y_label='Velocity (m/s)')

fdsplotlib.plot_to_fig(x_data=HRR_crit_exp_sup, y_data=V_crit_exp_sup, marker_style='kd', marker_fill_color='none', figure_handle=fig, data_label='Backlayering Controlled (Exp)')
fdsplotlib.plot_to_fig(x_data=HRR_crit_mod_sup, y_data=V_crit_mod_sup, marker_style='rd', marker_fill_color='none', figure_handle=fig, data_label='Backlayering Controlled (FDS)')
fdsplotlib.plot_to_fig(x_data=HRR_crit_exp_nosup, y_data=V_crit_exp_nosup, marker_style='ks', marker_fill_color='black', figure_handle=fig, data_label='Backlayering Not Controlled (Exp)')
fdsplotlib.plot_to_fig(x_data=HRR_crit_mod_nosup, y_data=V_crit_mod_nosup, marker_style='rs', marker_fill_color='red', figure_handle=fig, data_label='Backlayering Not Controlled (FDS)')

ax = plt.gca()
xtick = [10, 50, 100]
ax.set_xticks(xtick)
ax.set_xticklabels(['10','50','100'])

plt.savefig(pltdir + 'Critical_Velocity.pdf', format='pdf')
plt.close()


# Process Cold Flow case

M = importdata(outdir + 'Cold_Flow_Series_1_cat_devc.csv', ',', 2)
E = importdata(expdir + 'Cold_Flow_Series_1.csv', ',', 2)

mod_time_index = np.zeros(15, dtype=int)
for j in range(15):
    mod_time_index[j] = int(np.interp(300 * (j + 1), M.data[:, 0], np.arange(len(M.data[:, 0]))))

fig = fdsplotlib.plot_to_fig(x_data=[3,3,9,15,15], y_data=[169.4,164.7,292.6,372.4,379.9], marker_style='k^',
                             data_label='Measured',
                             x_min=0, x_max=16, y_min=0, y_max=400,
                             revision_label=version_string,
                             x_label='Number of Fans',
                             y_label='Volume Flow (m$^3$/s)')
fdsplotlib.plot_to_fig(x_data=M.data[mod_time_index,0]/300, y_data=M.data[mod_time_index,1], 
                       marker_style='ko-', marker_fill_color='none', data_label='FDS', figure_handle=fig)

ax = plt.gca()
ax.set_xticks([1, 3, 5, 7, 9, 11, 13, 15])

plt.savefig(pltdir + 'Cold_Flow_Series_1_Volume_Flow.pdf', format='pdf')
plt.close()

