# C. Paul (Orig. author: McGrattan)
# Converted to Python
# Original: 11-14-2022 | nat_conv_hot_plate.m

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import fdsplotlib

# -------------------------------
# Setup and constants
# -------------------------------
plot_style = fdsplotlib.get_plot_style('fds')
fds_dir = os.path.normpath(os.path.join(os.path.dirname(__file__),'..','..','..'))
results_dir = os.path.normpath(os.path.join(fds_dir,'..','out/Convection'))
pltdir = os.path.join(fds_dir, 'Manuals','FDS_Verification_Guide','SCRIPT_FIGURES','')
g = 9.80665
S = np.array([1.0, 0.1, 5.0, 10.0, 0.05])
T1 = np.array([503, 503, 503, 323, 323])
T2 = 293.0
MW = 28.85476       # FDS 'LJ AIR'
P0 = 101325.0
mu = 1.8216e-5
cp = 1000.0
k = 0.018216        # Pr=1 fluid

# -------------------------------
# Correlation curves
# -------------------------------
RAYLEIGH_DOWN = np.logspace(4, 9, 100)
RAYLEIGH_UP   = np.logspace(4, 11, 100)

NUSSELT_DOWN = 0.52 * RAYLEIGH_DOWN**0.2  # Eq. 9.32 (Incropera & DeWitt)
NUSSELT_UP = np.where(
    RAYLEIGH_UP < 1e7,
    0.54 * RAYLEIGH_UP**0.25,             # Eq. 9.30
    0.15 * RAYLEIGH_UP**0.333             # Eq. 9.31
)

# -------------------------------
# Plot setup
# -------------------------------
curves = [(RAYLEIGH_DOWN, NUSSELT_DOWN, 'k-', 'Correlation; Down')]
curves.append((RAYLEIGH_UP, NUSSELT_UP, 'k--', 'Correlation; Up'))
down_curves = []
up_curves = []

# -------------------------------
# FDS results
# -------------------------------
casename = [
    'nat_conv_hot_plate_1',
    'nat_conv_hot_plate_2',
    'nat_conv_hot_plate_3',
    'nat_conv_hot_plate_4',
    'nat_conv_hot_plate_5'
]
marker_style = ['gs', 'rs', 'ks', 'go', 'ro', 'ko']
res = ['8', '16', '32']

for j, r in enumerate(res):
    for i, s in enumerate(S):
        csv_file = os.path.join(results_dir, f'{casename[i]}_{r}_devc.csv')
        try:
            M = pd.read_csv(csv_file, header=1)
        except FileNotFoundError:
            print(f" Missing file: {csv_file}")
            continue
        Qdot_down = 1000.0 * M.iloc[-1, 1]
        Qdot_up   = 1000.0 * M.iloc[-1, 2]
        A = s**2
        Tm = 0.5 * (T1[i] + T2)
        beta = 1 / Tm
        rho = P0 * MW / (8341.5 * Tm)
        alpha = k / (rho * cp)
        nu = mu / rho
        L = s / 4
        Ra_FDS = (g * beta * (T1[i] - T2) * L**3) / (alpha * nu)
        Nu_FDS_down = (Qdot_down / A) * (L / k) / (T1[i] - T2)
        Nu_FDS_up   = (Qdot_up / A)   * (L / k) / (T1[i] - T2)
        if i==0:
            downLabel=f'Down ($S/\\delta x={r}$)'
            upLabel=f'Up ($S/\\delta x={r}$)'
        else:
            downLabel=""
            upLabel=""
        down_curves.append((Ra_FDS, Nu_FDS_down, marker_style[j], downLabel))
        up_curves.append((Ra_FDS, Nu_FDS_up, marker_style[j+3], upLabel))


curves.extend(down_curves)
curves.extend(up_curves)

# ============Plotting curves==================
xmax = max(np.max(x_data) for x_data, _, _, _ in curves)
ymax = max(np.max(y_data) for _, y_data, _, _ in curves)

git_file = os.path.join(results_dir, 'nat_conv_hot_plate_1_8_git.txt') 
version_string = fdsplotlib.get_version_string(git_file) if os.path.exists(git_file) else '' 

fig = fdsplotlib.plot_to_fig(x_data=RAYLEIGH_DOWN, y_data=np.ones_like(RAYLEIGH_DOWN), marker_style='k-',
                             markersize=11,revision_label=version_string,legend_location='upper left', 
                             x_label='Rayleigh Number', y_label='Nusselt Number',x_min=1, x_max=xmax, 
                             y_min=1, y_max=ymax,plot_type='loglog',xticks=[1e0, 1e5, 1e10, 1e15], 
                             yticks=[1, 2, 5, 10, 20, 50, 100, 200, 500, 1000]) 

for x_data, y_data, style, label in curves:
    fdsplotlib.plot_to_fig(x_data=x_data, y_data=y_data, figure_handle=fig, marker_style=style, data_label=label)

# ax = plt.gca() 
# ax.set_xticks([1e0, 1e5, 1e10, 1e15])
# ax.set_xticklabels([r'$10^0$', r'$10^5$', r'$10^{10}$', r'$10^{15}$'])
# ax.set_yticks([1, 2, 5, 10, 20, 50, 100, 200, 500, 1000])
# ax.set_yticklabels(['$1$', '$2$', '$5$', '$10$','$20$','$50$','$100$','$200$','$500$','$1000$'])

# -------------------------------
# Save figure
# -------------------------------
os.makedirs(pltdir, exist_ok=True)
out_pdf = os.path.join(pltdir, "nat_conv_hot_plate.pdf")
plt.savefig(out_pdf, format='pdf')
plt.close(fig)
