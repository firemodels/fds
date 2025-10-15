import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import fdsplotlib

plot_style = fdsplotlib.get_plot_style('fds')

# Paths and identifiers
outdir = '../../Verification/WUI/'
pltdir = '../../Manuals/FDS_Verification_Guide/SCRIPT_FIGURES/'

# Simulation parameters
U_vals = [0, 5]  # wind speeds (m/s)
R0 = 0.05  # zero-wind, zero-slope spread rate (m/s)
theta_w = 90  # wind compass direction (deg)
theta_s = 45  # compass direction of increasing slope (deg)
slope_angles = [0, 30]  # degrees of slope tested
theta_dev = np.arange(0, 360, 45)  # device angles

# Global parameters of fuel bed
sigma = 50.0 # svr (1/cm)
beta = 0.01 # packing ratio (-)

def wind_phi(U: float, theta_w_deg: float) -> tuple[float, float]:
    U_min = U * 60.0
    beta_op = 0.20395 * (sigma ** -0.8189)
    beta_ratio = beta / beta_op
    C = 7.47 * np.exp(-0.8711 * sigma ** 0.55)
    B = 0.15988 * sigma ** 0.54
    E = 0.715 * np.exp(-0.01094 * sigma)
    phi_w_x = C * (3.281 * U_min) ** B * (beta_ratio) ** -E * np.sin(np.deg2rad(theta_w_deg))
    phi_w_y = C * (3.281 * U_min) ** B * (beta_ratio) ** -E * np.cos(np.deg2rad(theta_w_deg))
    return phi_w_x, phi_w_y


def slope_phi(dzdx: float, dzdy: float) -> tuple[float, float]:
    dzds = np.sqrt(dzdx ** 2 + dzdy ** 2)
    coeff = 5.275 * beta ** -0.3
    phi_s_x = coeff * dzdx * dzds
    phi_s_y = coeff * dzdy * dzds
    return phi_s_y, phi_s_x


def virtual_wind(phi_s_x: float, phi_s_y: float) -> tuple[float, float]:
    beta_op = 0.20395 * (sigma ** -0.8189)
    beta_ratio = beta / beta_op
    C = 7.47 * np.exp(-0.8711 * sigma ** 0.55)
    B = 0.15988 * sigma ** 0.54
    E = 0.715 * np.exp(-0.01094 * sigma)
    phi_s = np.sqrt(phi_s_x ** 2 + phi_s_y ** 2)
    if phi_s > 0:
        uv_tmp = 0.3048 / phi_s * (phi_s / C * beta_ratio ** E) ** (1 / B)
    else:
        uv_tmp = 0.0
    u_virtual = uv_tmp * phi_s_x / 60.0
    v_virtual = uv_tmp * phi_s_y / 60.0
    return u_virtual, v_virtual


colors = ['k', 'r', 'b', 'g']
error_table = None

for slope in slope_angles:
    fig = None

    for idx_u, U in enumerate(U_vals):
        
        # Read arrival times
        CHID = f'LS_ellipse_{U}ms_{slope:02d}deg'
        devc_file = os.path.join(outdir, f'{CHID}_devc.csv')
        if not os.path.exists(devc_file):
            print(f'Error: File {devc_file} does not exist. Skipping case.')
            sys.exit()

        DEVC = pd.read_csv(devc_file, skiprows=1)

        # Distance from ignition for each point (along sloped surface)
        dist = np.sqrt(4.0 ** 2 + (4.0 * np.tan(np.deg2rad(slope)) * np.cos(np.deg2rad(theta_dev - theta_s))) ** 2)

        # Spread rate from arrival time at each point (use last row of DEVC)
        arrival_times = DEVC.iloc[-1, 1:].to_numpy(dtype=float)
        r = dist / arrival_times

        # Plot distance traveled in 10 seconds (Cartesian)
        x_pts = 10.0 * r * np.sin(np.deg2rad(theta_dev))
        y_pts = 10.0 * r * np.cos(np.deg2rad(theta_dev))

        # Revision string
        git_file = os.path.join(outdir, f'{CHID}_git.txt')
        version_string = fdsplotlib.get_version_string(git_file)
        
        if fig is None:
            fig = fdsplotlib.plot_to_fig(
                x_data=x_pts,
                y_data=y_pts,
                marker_style=f'{colors[idx_u]}o',
                x_min=-2,
                x_max=10,
                y_min=-4,
                y_max=5,
                revision_label=version_string,
                plot_title=f'Level set spread rate (LS_ellipse_*ms_{slope:02d}deg)',
                x_label='x (m)',
                y_label='y (m)',
                data_label=f'FDS {U} m/s'
            )
        else:
            fdsplotlib.plot_to_fig(x_data=x_pts,
                                   y_data=y_pts,
                                   marker_style=f'{colors[idx_u]}o',
                                   data_label=f'FDS {U} m/s',
                                   figure_handle=fig)

        # Theoretical spread
        HT = 0.2 # fuel layer height (m)
        UMF = U * 1.83 / np.log((20 + 1.18 * HT) / (0.43 * HT))

        phi_w_x, phi_w_y = wind_phi(UMF, theta_w)
        phi_s_x, phi_s_y = slope_phi(np.sin(np.deg2rad(theta_s)) * np.tan(np.deg2rad(slope)),
                                      np.cos(np.deg2rad(theta_s)) * np.tan(np.deg2rad(slope)))
        
        # combine slope and wind factors for total head spread rate
        phi_x = phi_w_x + phi_s_x
        phi_y = phi_w_y + phi_s_y
        R = R0 * (1.0 + np.sqrt(phi_x ** 2 + phi_y ** 2))

        # get virtual wind vector and apply to elliptical model
        u_v, v_v = virtual_wind(phi_s_x, phi_s_y)
        u_v = UMF * np.sin(np.deg2rad(theta_w)) + u_v
        v_v = UMF * np.cos(np.deg2rad(theta_w)) + v_v
        UMF_v = np.sqrt(u_v ** 2 + v_v ** 2)
        theta_w_v = 90.0 - np.degrees(np.arctan2(v_v, u_v))
        

        LB = 0.936 * np.exp(0.2566 * UMF_v) + 0.461 * np.exp(-0.1548 * UMF_v) - 0.397
        LB = float(np.clip(LB, 1.0, 8.0))
        HB = (LB + np.sqrt(LB ** 2 - 1.0)) / (LB - np.sqrt(LB ** 2 - 1.0))
        b = 0.5 * (R + R / HB)
        a = b / LB
        c = b - R / HB

        ang = np.arange(0, 360)
        x_s = np.cos(np.deg2rad(ang))
        y_s = np.sin(np.deg2rad(ang))

        twv = np.deg2rad(theta_w_v)
        D = (a ** 2 * (x_s * np.sin(twv) + y_s * np.cos(twv)) ** 2 +
             b ** 2 * (x_s * np.cos(twv) - y_s * np.sin(twv)) ** 2) ** -0.5
        R_x = D * (a ** 2 * np.cos(twv) * (x_s * np.sin(twv) + y_s * np.cos(twv)) -
                   b ** 2 * np.sin(twv) * (x_s * np.cos(twv) - y_s * np.sin(twv))) + c * np.sin(twv)
        R_y = D * (-a ** 2 * np.sin(twv) * (x_s * np.sin(twv) + y_s * np.cos(twv)) -
                   b ** 2 * np.cos(twv) * (x_s * np.cos(twv) - y_s * np.sin(twv))) + c * np.cos(twv)
        
        # theoretical spread distance after 10 s
        x_th = 10.0 * R_x
        y_th = 10.0 * R_y

        fdsplotlib.plot_to_fig(x_data=x_th,
                               y_data=y_th,
                               marker_style=f'{colors[idx_u]}-',
                               data_label=f'Expected {U} m/s',
                               figure_handle=fig)

        # Error computation for device angles
        R_target = np.sqrt(R_x ** 2 + R_y ** 2)
        theta_target = 90.0 - np.degrees(np.arctan2(R_y, R_x))
        theta_target = np.where(theta_target < 0, 360.0 + theta_target, theta_target)
        sort_idx = np.argsort(theta_target)
        theta_sorted = theta_target[sort_idx]
        R_sorted = R_target[sort_idx]
        R_interp = np.interp(theta_dev, theta_sorted, R_sorted, left=None, right=None)
        err_vals = np.abs(r - R_interp) / R_interp

        col_name = fr'{U} m/s; {slope:d}$^\circ$'
        if error_table is None:
            error_table = pd.DataFrame({'angle': theta_dev, col_name: err_vals})
        else:
            error_table[col_name] = err_vals

    # Plot ignition point
    ax = plt.gca()
    ax.plot(0.0, 0.0, 'rx', label='_nolegend_')

    # Add wind annotation arrow
    ax.annotate('wind', xytext=(2.0, -3.5), xy=(4.0, -3.5), xycoords='data',
                ha='center', va='center',
                arrowprops=dict(arrowstyle='->', color='k', lw=1.0)) 
    # Add slope annotation arrow (only if slope > 0)
    if slope > 0:
        dx = 2.0 * np.cos(np.deg2rad(theta_s))
        dy = 2.0 * np.sin(np.deg2rad(theta_s))
        ax.annotate('slope', xytext=(2.0, -2.75), xy=(2.0 + dx, -2.75 + dy),
                    ha='center', va='center',
                    arrowprops=dict(arrowstyle='->', color='k', lw=1.0))

    # Save per-slope figure
    fig.savefig(os.path.join(pltdir, f'level_set_ellipse_{slope:02d}deg.pdf'), format='pdf')
    plt.close(fig)


# Tolerance check
if error_table is not None:
    max_err = float(np.max(error_table.drop(columns=['angle']).to_numpy()))
    if max_err > 0.2:
        print(f'Python Warning: LS_ellipse is out of tolerance. Max error = {max_err:.3f}')

    # Error table plot
    cols = [c for c in error_table.columns if c != 'angle']
    # Initialize figure with first series
    fig_err = fdsplotlib.plot_to_fig(x_data=error_table['angle'].to_numpy(),
                                     y_data=error_table[cols[0]].to_numpy(),
                                     marker_style=f'{colors[0]}-',
                                     revision_label=version_string,
                                     plot_title='Level set spread error',
                                     x_label=r'compass angle ($^\circ$)',
                                     y_label='relative error (-)',
                                     y_min=0,
                                     y_max=0.3,
                                     data_label=cols[0])
    
    for i, c in enumerate(cols[1:], start=1):
        fdsplotlib.plot_to_fig(x_data=error_table['angle'].to_numpy(),
                               y_data=error_table[c].to_numpy(),
                               marker_style=f'{colors[i % len(colors)]}-',
                               data_label=c,
                               figure_handle=fig_err)



    fig_err.savefig(os.path.join(pltdir, 'level_set_ellipse_error.pdf'), format='pdf')
    plt.close(fig_err)


