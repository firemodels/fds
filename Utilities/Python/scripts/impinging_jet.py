#!/usr/bin/env python3
"""
impinging_jet.py
Converted from MATLAB script impinging_jet.m to Python.

Generates:
  - impinging_jet_correlation.pdf
  - impinging_jet_local_Re_1e5.pdf
  - impinging_jet_local_Re_4e5.pdf
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
import fdsplotlib

def set_ticks_like_matlab(fig):
    ax = fig.axes[0]
    ax.tick_params(axis="both", direction="in", top=True, right=True, width=0.5)

    for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_linewidth(0.5)
    yticks = ax.yaxis.get_major_ticks()
    yticks[0].tick1line.set_visible(False)
    yticks[0].tick2line.set_visible(False)
    yticks[-1].tick1line.set_visible(False)
    yticks[-1].tick2line.set_visible(False)
    
    xticks = ax.xaxis.get_major_ticks()
    xticks[0].tick1line.set_visible(False)
    xticks[0].tick2line.set_visible(False)
    xticks[-1].tick1line.set_visible(False)
    xticks[-1].tick2line.set_visible(False)

plot_style = fdsplotlib.get_plot_style("fds")

# -----------------------------
# Paths and identifiers
# -----------------------------
firemodels_dir = os.path.normpath(os.path.join(os.path.dirname(__file__),'..','..','..','..'))
fds_dir = os.path.normpath(os.path.join(os.path.dirname(__file__),'..','..','..'))
outdir = os.path.join(firemodels_dir, 'out','Convection','')
pltdir = os.path.join(fds_dir, 'Manuals','FDS_Verification_Guide','SCRIPT_FIGURES','')
os.makedirs(pltdir, exist_ok=True)

chid = 'impinging_jet'

# -----------------------------
# Constants / inputs
# -----------------------------
D_h   = 0.2          # hydraulic diameter of the jet [m]
H     = 1.0          # distance from jet to wall [m]
mu    = 1.822e-05    # dynamic viscosity [kg/m/s]
k     = 2.566e-02    # thermal conductivity [W/m/K]
Pr    = 0.71         # Prandtl number [-]
T_j   = 100.0        # temperature of jet exit [C]
T_w   = 20.0         # constant plate temperature [C]
rho_a = 1.200        # ambient density [kg/m^3]
rho_j = rho_a * (T_w + 273.0) / (T_j + 273.0)  # jet fluid density [kg/m^3]
nu    = mu / rho_j   # kinematic viscosity [m^2/s]
U_j   = [10.0, 40.0] # jet exit velocities [m/s] corresponding to Re ~ 1e5 and 4e5
A_r   = 0.01         # D_h^2/(4 r^2), where r defines outer extent of the averaging region

res_list = ['coarse', 'medium', 'fine']
Re_str   = ['1e5', '4e5']

# Relative error tolerance for Nu comparison
E_tol = 0.1

# Get git hash
git_file = os.path.join(outdir, 'impinging_jet_Re_1e5_coarse_git.txt')
version_string = fdsplotlib.get_version_string(git_file)

# -----------------------------
# Plot 1: Martin correlation vs Re and FDS points
# -----------------------------
# Martin correlation coefficient G for given geometry
G = 2.0 * np.sqrt(A_r) * ((1 - 2.2*np.sqrt(A_r)) / (1 + 0.2*(H/D_h - 6)*np.sqrt(A_r)))

RE = np.linspace(2e3, 5e5, 1000)
NU = G * 2*np.sqrt(RE) * np.sqrt(1 + 0.005*RE**0.55) * (Pr**0.42)

fig = fdsplotlib.plot_to_fig(x_data=RE, y_data=NU, marker_style='k-', data_label='Martin', linewidth=2,
                             x_min=0.0, x_max=5e5, y_min=0, y_max=600,
                             revision_label=version_string,
                             x_label='Re',
                             y_label='Nu')


markers = {'coarse': ('b', 's'), 'medium': ('r', 's'), 'fine': ('g', 's')}

E_FDS = np.zeros((2, 3))  # shape: (len(Re_str), len(res_list))

# Loop over resolutions and Re cases
for j, res in enumerate(res_list):
    for i, Re_tag in enumerate(Re_str):
        # Compute Re from U_j for this case
        Re_j = D_h * U_j[i] / nu

        # Martin Nu for this Re
        Nu_cor = G * 2*np.sqrt(Re_j) * np.sqrt(1 + 0.005*Re_j**0.55) * (Pr**0.42)

        # FDS results: read *_devc.csv and get average HF over second half
        devc_path = os.path.join(outdir, f'impinging_jet_Re_{Re_tag}_{res}_devc.csv')
        df = pd.read_csv(devc_path, header=1)
        
        # Mean HF over second half of samples
        n = len(df)
        hf_col = 'HF' if 'HF' in df.columns else df.columns[-1]  # simple fallback
        HF_mean = df[hf_col].iloc[n // 2 :].mean()

        qconv = HF_mean * 1000.0  # W/m^2
        h_fds = qconv / (T_j - T_w)
        Nu_fds = h_fds * D_h / k

        # relative error
        E = abs(Nu_cor - Nu_fds) / abs(Nu_cor)
        E_FDS[i, j] = E
        if E > E_tol:
            print(f'MATLAB Warning (ported): impinging jet error = {E:.3f} at Re_j={Re_tag}, Res={res}')

        label = f'FDS {res}' if i == 0 else None
        fig = fdsplotlib.plot_to_fig(x_data=Re_j, y_data=Nu_fds, data_label=label, figure_handle=fig, marker_style=markers[res], marker_fill_color=(1,1,1,0.0), markeredgewidth=2.0, markersize=4)

# Add legend and remove the line from the markers added by fdsplotlib
fig.axes[0].legend(loc='upper left', frameon=True, fontsize=plot_style['Key_Font_Size'])
for handle in fig.gca().get_legend().legend_handles[1:]:
    handle.set_linestyle('')  # remove line
    
formatter = ScalarFormatter(useMathText=True)
formatter.set_powerlimits((0, 0))  # force scientific notation
fig.axes[0].xaxis.set_major_formatter(formatter)
fig.axes[0].xaxis.get_offset_text().set_fontsize(plot_style['Label_Font_Size'])
set_ticks_like_matlab(fig)

# Save figure
corr_pdf = os.path.join(pltdir, 'impinging_jet_correlation.pdf')
fig.savefig(corr_pdf)
plt.close(fig)

# -----------------------------
# Plot 2: Local Nu(x) profiles for each Re case
# -----------------------------
styles = ['k-.', 'k--', 'k-']  # coarse, medium, fine
dashes = [(12, 7.4, 1, 7.4), (18.0, 11.1), (None, None)]
for i, Re_tag in enumerate(Re_str):
    fig = None
    for j, res in enumerate(res_list):
        line_path = os.path.join(outdir, f'impinging_jet_Re_{Re_tag}_{res}_line.csv')
        LINE = pd.read_csv(line_path, header=1)

        # x is first column; QCONV column used (match MATLAB)
        x = LINE.iloc[:, 0].to_numpy()
        q_col = 'QCONV' if 'QCONV' in LINE.columns else LINE.columns[-1]
        q_x = LINE[q_col].to_numpy() * 1000.0  # W/m^2

        Nu_x = (q_x / (T_j - T_w)) * D_h / k

        revlabel = version_string if j == 0 else None
        fig = fdsplotlib.plot_to_fig(x_data=x, y_data=Nu_x, marker_style=styles[j], data_label=r'$D_{h}/\delta$x = %d'%(7*(2**j)), linewidth=0.5,
                                     x_min=-0.5, x_max=0.5, y_min=0, y_max=1000,
                                     revision_label=revlabel, figure_handle=fig, 
                                     x_label='x (m)',
                                     y_label=r'Nu$_{D_h}$(x)', line_dashes=dashes[j])
        
    fig.axes[0].legend(loc='upper left', frameon=True, fontsize=plot_style['Key_Font_Size'])
    set_ticks_like_matlab(fig)
    fig.axes[0].set_xticks([-0.5,0,0.5])
    fig.axes[0].xaxis.set_tick_params(pad=15) 
    local_pdf = os.path.join(pltdir, f'impinging_jet_local_Re_{Re_tag}.pdf')
    fig.savefig(local_pdf)
    plt.close(fig)


print('Completed plots:')
print(' -', corr_pdf)
for Re_tag in Re_str:
    print(' -', os.path.join(pltdir, f'impinging_jet_local_Re_{Re_tag}.pdf'))
