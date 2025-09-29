
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# include FDS plot styles, etc.
import fdsplotlib

# Get plot style parameters
plot_style = fdsplotlib.get_plot_style('fds')

expdir = '../../../exp/Submodules/macfp-db/Gaseous_Pool_Fires/McCaffrey_Flames/Computational_Results/2017/Data/'
outdir = '../../../out/'
Manuals_Dir = '../../Manuals/'
datadir = os.path.join(outdir, 'McCaffrey_Plume')
pltdir = os.path.join(Manuals_Dir, 'FDS_Validation_Guide', 'SCRIPT_FIGURES', 'McCaffrey_Plume')
git_file =  os.path.join(datadir, f'McCaffrey_14_kW_11_git.txt')
version_string = fdsplotlib.get_version_string(git_file)

Q = np.array([14.4, 21.7, 33.0, 44.9, 57.5])  # kW
chid = ['McCaffrey_14kW', 'McCaffrey_22kW', 'McCaffrey_33kW', 'McCaffrey_45kW', 'McCaffrey_57kW']
mark = ['ko', 'k+', 'k^', 'ks', 'kd']

# Velocity and temperature correlations

zq = np.logspace(-2, 0, 100)
vq = np.zeros_like(zq)
Tq = np.zeros_like(zq)
for i, z in enumerate(zq):
    if z < 0.08:
        vq[i] = 6.84 * z**0.5
        Tq[i] = 800.0
    elif 0.08 <= z <= 0.2:
        vq[i] = 1.93
        Tq[i] = 63.0 * z**(-1)
    else:
        vq[i] = 1.12 * z**(-1.0/3.0)
        Tq[i] = 21.6 * z**(-5.0/3.0)

# Plot velocity correlation and measurement points

fig1 = fdsplotlib.plot_to_fig(x_data=zq, y_data=vq, data_label='$(z/Q^{2/5})^\\eta$', 
                              marker_style='b--',
                              x_min=0.01, x_max=1, y_min=0.5, y_max=3,
                              x_label='$z/Q^{2/5}$ $(\\hbox{m kW}^{-2/5})$',
                              y_label='${V/Q}^{1/5}$ (m s$^{-1}$ kW$^{-1/5}$)',
                              legend_location='lower right',
                              plot_type='loglog')

for i, c in enumerate(chid):
    fname = os.path.join(expdir, f'{c}_V.csv')
    if os.path.exists(fname):
        M = pd.read_csv(fname, comment=None)
        fdsplotlib.plot_to_fig(x_data=M.iloc[:,0].values,y_data=M.iloc[:, 1].values, 
                               figure_handle=fig1,
                               marker_style=mark[i], data_label=f'{Q[i]:.1f} kW')

ax = plt.gca()
ax.text(0.012,1.2, '$\\eta=1/2$' , fontsize=plot_style['Label_Font_Size'])
ax.text(0.100,2.3, '$\\eta=0$'   , fontsize=plot_style['Label_Font_Size'])
ax.text(0.400,1.8, '$\\eta=-1/3$', fontsize=plot_style['Label_Font_Size'])

fig1.savefig(os.path.join(pltdir, 'McCaffrey_Velocity_Correlation.pdf'), format='pdf')
plt.close(fig1)

# Plot temperature correlation and measurement points

fig2 = fdsplotlib.plot_to_fig(x_data=zq, y_data=Tq, data_label='$(z/Q^{2/5})^\\eta$',
                              marker_style='r--',
                              x_min=0.008, x_max=1, y_min=100, y_max=1200,
                              x_label='$z/Q^{2/5}$ $(\\hbox{m kW}^{-2/5})$',
                              y_label='$\\Delta T \\; (^\\circ$C)',
                              legend_location='lower left',
                              plot_type='loglog')

for i, c in enumerate(chid):
    fname = os.path.join(expdir, f'{c}_T.csv')
    if os.path.exists(fname):
        M = pd.read_csv(fname, comment=None)
        fdsplotlib.plot_to_fig(x_data=M.iloc[:,0].values,y_data=M.iloc[:, 1].values,
                               figure_handle=fig2,
                               marker_style=mark[i], data_label=f'{Q[i]:.1f} kW')

ax = plt.gca()
ax.text(0.012,500, '$\\eta=0$'   , fontsize=plot_style['Label_Font_Size'])
ax.text(0.060,400, '$\\eta=-1$'  , fontsize=plot_style['Label_Font_Size'])
ax.text(0.100,150, '$\\eta=-5/3$', fontsize=plot_style['Label_Font_Size'])

fig2.savefig(os.path.join(pltdir, 'McCaffrey_Temperature_Correlation.pdf'), format='pdf')
plt.close(fig2)

# Extrapolate temperature correlation to burner surface

os.makedirs(pltdir, exist_ok=True)
for i, c in enumerate(chid):
    fname = os.path.join(expdir, f'{c}_T.csv')
    if not os.path.exists(fname):
        continue
    M = pd.read_csv(fname, comment=None)
    ZS = M.iloc[:, 0].values
    TC = M.iloc[:, 1].values

    fig = fdsplotlib.plot_to_fig(x_data=zq, y_data=Tq, data_label='$(z/Q^{2/5})^\\eta$',
                                 marker_style='r--',
                                 x_min=0.008, x_max=1, y_min=100, y_max=1200,
                                 x_label='$z/Q^{2/5}$ $(\\hbox{m kW}^{-2/5})$',
                                 y_label='$\\Delta T \\; (^\\circ$C)',
                                 legend_location='lower left',
                                 plot_title=f'McCaffrey Centerline Temperature Data, {Q[i]:.1f} kW',
                                 plot_type='loglog')

    fdsplotlib.plot_to_fig(x_data=ZS,y_data=TC,
                           figure_handle=fig,
                           marker_style=mark[i], data_label=f'{Q[i]:.1f} kW')

    # identify near-surface (ZS < 0.05)
    isurf = np.where(ZS < 0.05)[0]
    if len(isurf) >= 2:
        fdsplotlib.plot_to_fig(x_data=ZS[isurf], y_data=TC[isurf], marker_style='bo', figure_handle=fig)
        # least squares fit (linear in this small z region)
        A = np.column_stack((ZS[isurf], np.ones_like(ZS[isurf])))
        b = TC[isurf]
        x, *_ = np.linalg.lstsq(A, b, rcond=None)
        T_surf = x[1]
        z1 = 0.008; z2 = 0.05
        T1 = x[0]*z1 + x[1]
        T2 = x[0]*z2 + x[1]
        fdsplotlib.plot_to_fig(x_data=[z1,z2], y_data=[T1,T2], marker_style='g-', figure_handle=fig, data_label='Least Squares Fit')
        ax = plt.gca()
        ax.text(0.009, 400., '$T_{\\rm surf}$='+f'{T_surf:3.0f} $^\\circ$C', fontsize=plot_style['Label_Font_Size'])
    else:
        T_surf = np.nan

    fig.savefig(os.path.join(pltdir, f'{c}_Surface_Temp.pdf'), format='pdf')
    plt.close(fig)

# Plot FDS predictions

T0 = 273.15 + 20  # K

resolution = ['Crude', 'Coarse', 'Medium', 'Fine']
for res in resolution:
    if res == 'Crude':
        chid = ['McCaffrey_14_kW_5','McCaffrey_22_kW_5','McCaffrey_33_kW_5','McCaffrey_45_kW_5','McCaffrey_57_kW_5']
    elif res == 'Coarse':
        chid = ['McCaffrey_14_kW_11','McCaffrey_22_kW_11','McCaffrey_33_kW_11','McCaffrey_45_kW_11','McCaffrey_57_kW_11']
    elif res == 'Medium':
        chid = ['McCaffrey_14_kW_21','McCaffrey_22_kW_21','McCaffrey_33_kW_21','McCaffrey_45_kW_21','McCaffrey_57_kW_21']
    elif res == 'Fine':
        chid = ['McCaffrey_14_kW_45','McCaffrey_22_kW_45','McCaffrey_33_kW_45','McCaffrey_45_kW_45','McCaffrey_57_kW_45']

    # Velocity correlation figures (one for each resolution)

    fig_v = fdsplotlib.plot_to_fig(x_data=zq, y_data=vq, data_label='$(z/Q^{2/5})^\\eta$',
                                   marker_style='b--',
                                   x_min=0.01, x_max=1, y_min=0.5, y_max=3,
                                   revision_label=version_string,
                                   x_label='$z/Q^{2/5}$ $(\\hbox{m kW}^{-2/5})$',
                                   y_label='$V/Q^{1/5}$ $(\\hbox{m s}^{-1} \\; \\hbox{kW}^{-1/5})$',
                                   plot_title=f'McCaffrey Centerline Velocity, {res}',
                                   legend_location='lower right',
                                   plot_type='loglog')

    for i, c in enumerate(chid):
        fname = os.path.join(datadir, f'{c}_line.csv')
        if os.path.exists(fname):
            M = pd.read_csv(fname, header=1)
            # find columns by name if present; else use first/second
            if 'Height' in M.columns and 'vel' in M.columns:
                z = M['Height'].values
                v = M['vel'].values
            else:
                z = M.iloc[:, 0].values
                v = M.iloc[:, 1].values
            zq_fds = z / (Q[i]**(2.0/5.0))
            vq_fds = v / (Q[i]**(1.0/5.0))
            fdsplotlib.plot_to_fig(x_data=zq_fds, y_data=vq_fds, 
                                          data_label=f'{Q[i]} kW', marker_style=mark[i], figure_handle=fig_v)
    fig_v.savefig(os.path.join(pltdir, f'McCaffrey_Velocity_Correlation_{res}.pdf'))
    plt.close(fig_v)

    # Temperature correlation figures (one for each resolution)

    fig_t = fdsplotlib.plot_to_fig(x_data=zq, y_data=Tq, data_label='$(z/Q^{2/5})^\\eta$',
                                   marker_style='r--',
                                   x_min=0.008, x_max=1, y_min=100, y_max=1200,
                                   revision_label=version_string,
                                   x_label='$z/Q^{2/5}$ $(\\hbox{m kW}^{-2/5})$',
                                   y_label='$\\Delta{T}$ ($^\\circ$C)',
                                   plot_title=f'McCaffrey Centerline Temperature, {res}',
                                   legend_location='lower right',
                                   plot_type='loglog')

    for i, c in enumerate(chid):
        fname = os.path.join(datadir, f'{c}_line.csv')
        if os.path.exists(fname):
            M = pd.read_csv(fname, header=1)
            if 'Height' in M.columns and 'tmp' in M.columns:
                z = M['Height'].values
                T = M['tmp'].values + 273.15
            else:
                z = M.iloc[:, 0].values
                T = M.iloc[:, 2].values + 273.15 if M.shape[1] > 2 else M.iloc[:, 1].values + 273.15
            zq_fds = z / (Q[i]**(2.0/5.0))
            fdsplotlib.plot_to_fig(x_data=zq_fds, y_data=T-T0, 
                                          data_label=f'{Q[i]} kW', marker_style=mark[i], figure_handle=fig_t)
    fig_t.savefig(os.path.join(pltdir, f'McCaffrey_Temperature_Correlation_{res}.pdf'))
    plt.close(fig_t)

