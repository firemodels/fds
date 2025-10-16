#!/usr/bin/env python
"""
mass_balance.py
Validates mass conservation comparing time derivative of species mass with inlet/outlet flow rates.
Original MATLAB script by McDermott (05 Dec 2017)
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import fdsplotlib
import os

def plot_mass_balance(chid, title_text):
    
    plot_style = fdsplotlib.get_plot_style('fds')
    ddir   = '../../Verification/Species/'
    pltdir = '../../Manuals/FDS_Verification_Guide/SCRIPT_FIGURES/'
    
    mass_file = os.path.join(ddir, f'{chid}_mass.csv')

    try:
        M = pd.read_csv(mass_file, skiprows=1)  # Skip units row, use column names as header
    except FileNotFoundError:
        print(f'Error: File {mass_file} does not exist. Skipping case.')
        return 0

    t = M.iloc[:, 0].values
    m = M['WATER VAPOR'].values
    
    # Compute dm/dt using numerical differentiation
    dmdt = np.zeros(len(t))
    for i in range(1, len(t)):
        dmdt[i] = (m[i] - m[i-1]) / (t[i] - t[i-1])
    
    devc_file = os.path.join(ddir, f'{chid}_devc.csv')
    
    F = pd.read_csv(devc_file, skiprows=1)  # Skip units row, use column names as header
    mdot_in  = F['H2O in'].values
    mdot_out = F['H2O out'].values
    
    # Residual: dm/dt - mdot_in - mdot_out (should be ~0)
    bal = dmdt - mdot_in - mdot_out
        
    # Get version string if git file exists
    git_file = os.path.join(ddir, f'{chid}_git.txt')
    version_string = fdsplotlib.get_version_string(git_file) if os.path.exists(git_file) else ''
    
    # Plot zero reference line first (black solid line, no label)
    fig = fdsplotlib.plot_to_fig(x_data=t, y_data=np.zeros_like(t), marker_style='k-',revision_label=version_string,legend_location='upper right',
                                 x_label='time (s)', y_label='mass flow (kg/s)',x_min=0, x_max=2000, y_min=-0.005, y_max=0.015)
    fdsplotlib.plot_to_fig(x_data=t, y_data=mdot_in, figure_handle=fig, marker_style='r-', data_label='Inlet H2O')
    fdsplotlib.plot_to_fig(x_data=t, y_data=-mdot_out, figure_handle=fig, marker_style='m-', data_label='Outlet H2O')
    fdsplotlib.plot_to_fig(x_data=t, y_data=bal, figure_handle=fig, marker_style='b-', data_label='dm/dt+out-in')
    
    ax    = plt.gca()
    lines = ax.get_lines()
    lines[2].set_color('#EDB120')         # Orange for outlet (second data line)
    ax.locator_params(axis='x', nbins=6)  # ~6 ticks on x-axis
    ax.locator_params(axis='y', nbins=5)  # ~5 ticks on y-axis
    fdsplotlib.apply_global_exponent(ax, axis='y')   
    ax.text(100, 18e-3, title_text,fontsize=plot_style['Title_Font_Size'],fontname=plot_style['Font_Name'])
    
    # Compute mean error for t > 1000 s
    mask = t > 1000
    if np.any(mask):
        mass_error = np.abs(np.mean(bal[mask]))
        if mass_error > 1e-5:
            print(f'Python Warning: mass error = {mass_error:.6e} in {chid}')
    else:
        mass_error = np.nan
        print(f'Warning: No data points with t > 1000 s in {chid}')
    
    # Save figure
    os.makedirs(pltdir, exist_ok=True)
    output_file = os.path.join(pltdir, f'{chid}_mass_balance.pdf')
    plt.savefig(output_file, format='pdf')
    plt.close(fig)
        
    return mass_error

if __name__ == '__main__':
    """Main execution block"""
    
    print("Running mass balance verification plots...")
    error1 = plot_mass_balance('mass_flux_wall_yindex','Primitive Species Mass Balance')
    error2 = plot_mass_balance('mass_flux_wall_zindex','Lumped Species Mass Balance')
    print("Mass balance verification completed successfully!")