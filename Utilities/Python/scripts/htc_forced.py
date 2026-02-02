# C. Paul (Orig. author: McDermott)
# Converted to Python
# Original: 11-14-2022 | htc_forced.m

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import fdsplotlib

# === Directories ===
plot_style = fdsplotlib.get_plot_style('fds')
fds_dir = os.path.normpath(os.path.join(os.path.dirname(__file__),'..','..','..'))
outdir = os.path.normpath(os.path.join(fds_dir,'..','out/Convection'))
pltdir = os.path.join(fds_dir, 'Manuals','FDS_Verification_Guide','SCRIPT_FIGURES','')
chid = 'forced_conv_flat_plate'



# === Velocity and resolution configurations ===
vel = ['u0p1', 'u1', 'u10']
res = ['25cm', '10cm', '2p5cm']
res_style = ['b*', 'ks', 'mo']
u = [0.1, 1, 10]
dx = [0.25, 0.1, 0.025]
dx_style = ['b--', 'k--', 'm--']

# === Physical parameters ===
N = 32
L = 16
rho = 1.165  # kg/m^3
mu = 1e-5    # Pa·s
k = 1e-2     # W/m·K
T_w = 30     # °C
T_g = 20     # °C

x = np.linspace(0, L, N)
REL_ERROR = np.zeros((len(vel), len(res)))




# === Main Loop ===
for i in range(len(vel)):
    Re_x = (rho / mu) * u[i] * x
    Nu_x = 0.0296 * Re_x ** 0.8

    # Store all curves in a list of tuples: (x_data, y_data, style, label)
    curves = [(Re_x, Nu_x, 'k-', 'forced local')]
    fds_nu_curves = []
    fds_ref_curves = []

    for j in range(len(res)):
        # === Load FDS CSV ===
        fname = os.path.join(outdir, f"{chid}_{vel[i]}_{res[j]}_line.csv")
        try:
            M = pd.read_csv(fname, header=1)
        except FileNotFoundError:
            print(f" Missing file: {fname}")
            continue
        if not all(c in M.columns for c in ['HTC-x', 'HTC', 'QCONV']):
            print(f" Skipping invalid file: {fname}")
            continue

        x_fds = M['HTC-x'].values
        h_fds = M['HTC'].values
        q_fds = M['QCONV'].values

        DTMP_fds = -1000 * q_fds / h_fds
        Nu_x_fds = h_fds * DTMP_fds / (T_w - T_g) * x_fds / k
        Re_x_fds = (rho / mu) * u[i] * x_fds

        # Append FDS curves
        fds_nu_curves.append((Re_x_fds, Nu_x_fds, res_style[j], fr'FDS $\delta x$={res[j]}'))
        fds_ref_curves.append((Re_x_fds, 2 / dx[j] * x_fds, dx_style[j], fr'$2*1/{dx[j]}*x$'))

        # Compute relative error
        REL_ERROR[i, j] = abs(Nu_x[-1] - Nu_x_fds[-1]) / Nu_x[-1]


    curves.extend(fds_nu_curves)
    curves.extend(fds_ref_curves)
    
    # ============Plotting curves==================
    xmax = max(np.max(x_data) for x_data, _, _, _ in curves)
    ymax = max(np.max(y_data) for _, y_data, _, _ in curves)
    
    git_file = os.path.join(outdir, f'{chid}_{vel[i]}_{res[0]}_git.txt') 
    version_string = fdsplotlib.get_version_string(git_file) if os.path.exists(git_file) else '' 
    
    fig = fdsplotlib.plot_to_fig(x_data=Re_x, y_data=np.zeros_like(Re_x), marker_style='k-',
                                 revision_label=version_string,legend_location='upper left', 
                                 x_label=r'Re$_x$', y_label=r'Nu$_x$',x_min=0, x_max=xmax, 
                                 y_min=0, y_max=ymax) 
    # End the x and y axis at a tick location
    ax = plt.gca() 

    ticks = ax.xaxis.get_majorticklocs()
    if len(ticks) > 1:
        xmax_aligned = ticks[-1] if ticks[-1] >= xmax else ticks[-1] + (ticks[-1] - ticks[-2])
    else:
        xmax_aligned = xmax
    
    ticks = ax.yaxis.get_majorticklocs()
    if len(ticks) > 1:
        ymax_aligned = ticks[-1] if ticks[-1] >= ymax else ticks[-1] + (ticks[-1] - ticks[-2])
    else:
        ymax_aligned = ymax
    plt.close(fig) 
    
    # Create another one with corect xmax_aligned and ymax_aligned
    fig = fdsplotlib.plot_to_fig(x_data=Re_x, y_data=np.zeros_like(Re_x), marker_style='k-',
                                 revision_label=version_string,legend_location='upper left', 
                                 x_label=r'Re$_x$', y_label=r'Nu$_x$',x_min=0, x_max=xmax_aligned, 
                                 y_min=0, y_max=ymax_aligned) 
    ax = plt.gca() 
    ax.ticklabel_format(style='sci', axis='x', scilimits=(-4, 4))
    ax.ticklabel_format(style='sci', axis='y', scilimits=(-4, 4))
    
    # Plot title
    plotTitle=f'Velocity = {u[i]} m/s' 
    ax.text(xmax_aligned*0.5, ymax_aligned*0.9, plotTitle,fontsize=plot_style['Title_Font_Size'],fontname=plot_style['Font_Name'])

    for x_data, y_data, style, label in curves:
        fdsplotlib.plot_to_fig(x_data=x_data, y_data=y_data, figure_handle=fig, marker_style=style, data_label=label)

    # === Save figure ===
    os.makedirs(pltdir, exist_ok=True)
    output_file = os.path.join(pltdir, f"{chid}_{vel[i]}.pdf")
    plt.savefig(output_file, format='pdf')
    plt.close(fig)
    
    
# === Error test ===
REF_ERROR = np.array([
    [0.13117575, 0.16006822, 0.0400888],
    [0.30421425, 0.07066064, 0.24389138],
    [0.41742345, 0.22269597, 0.07470432]
])

# print(REL_ERROR)

ERROR = np.linalg.norm(REL_ERROR - REF_ERROR)

if ERROR > 0.1:
    print(f"Warning: Forced convection out of tolerance. ERROR = {ERROR:.3f}")

