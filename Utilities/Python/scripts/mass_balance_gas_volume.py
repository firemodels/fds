"""
mass_balance_gas_volume.py
Validates mass conservation in a gas volume comparing time derivative of mass with net mass flux through boundaries.
Original MATLAB script by McDermott (2 Dec 2021)
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import fdsplotlib
import os

def plot_mass_balance(chid, title_text):
    
    plot_style = fdsplotlib.get_plot_style('fds')
    fds_dir    = os.path.normpath(os.path.join(os.path.dirname(__file__),'..','..','..'))
    ddir       = os.path.join(fds_dir, 'Verification','Species','')
    pltdir     = os.path.join(fds_dir, 'Manuals','FDS_Verification_Guide','SCRIPT_FIGURES','')
    devc_file  = os.path.join(ddir, f'{chid}_devc.csv')

    try:
        M = pd.read_csv(devc_file, skiprows=1)  # Skip units row, use column names as header
    except FileNotFoundError:
        print(f'Error: File {devc_file} does not exist. Skipping case.')
        return 0

    t = M.iloc[:, 0].values
    m = M['M'].values
    dmdt = np.zeros(len(t))
    for i in range(1, len(t)):
        dmdt[i] = (m[i] - m[i-1]) / (t[i] - t[i-1])
    
    mf_x1 = M['TMF_X1'].values
    mf_x2 = M['TMF_X2'].values
    mf_y1 = M['TMF_Y1'].values
    mf_y2 = M['TMF_Y2'].values
    mf_z1 = M['TMF_Z1'].values
    mf_z2 = M['TMF_Z2'].values    
    # Net flux: out - in
    net_flux = mf_x2 - mf_x1 + mf_y2 - mf_y1 + mf_z2 - mf_z1
    # Residual: dm/dt + (out - in) (should be ~0)
    bal = dmdt + net_flux
        
    # Get version string if git file exists
    git_file = os.path.join(ddir, f'{chid}_git.txt')
    version_string = fdsplotlib.get_version_string(git_file) if os.path.exists(git_file) else ''
    
    # Plot zero reference line first (black solid line, no label)
    fig = fdsplotlib.plot_to_fig(x_data=t, y_data=np.zeros_like(t), marker_style='k-',
                                 revision_label=version_string, legend_location='upper right',x_label='time (s)', y_label='mass flow (kg/s)',
                                 x_min=0, x_max=2000, y_min=-0.004, y_max=0.004,xnumticks=5, ynumticks=5)
    fdsplotlib.plot_to_fig(x_data=t, y_data=dmdt, figure_handle=fig, marker_style='r-', data_label='dm/dt')
    fdsplotlib.plot_to_fig(x_data=t, y_data=net_flux, figure_handle=fig, marker_style='m-', data_label='out-in')
    fdsplotlib.plot_to_fig(x_data=t, y_data=bal, figure_handle=fig, marker_style='b-', data_label='dm/dt+out-in')
    
    ax = plt.gca()
    lines = ax.get_lines()
    lines[2].set_color('orange')  # Orange for net flux (second data line)
    legend = ax.get_legend()
    if legend:
        legend.legend_handles[1].set_color('orange')  # Second legend entry (out-in)
    fdsplotlib.apply_global_exponent(ax, axis='y')
    if title_text:
        ax.text(100, 18e-3, title_text, fontsize=plot_style['Title_Font_Size'], 
                fontname=plot_style['Font_Name'])
    
    # Check balance and report error
    mass_error = np.max(np.abs(bal))
    if mass_error > 1e-6:
        print(f'Python Warning: mass error = {mass_error:.6e} in {chid}')
    
    # Save figure
    output_file = os.path.join(pltdir, f'{chid}.pdf')
    plt.savefig(output_file, format='pdf')
    plt.close(fig)
        
    return mass_error

if __name__ == '__main__':
    """Main execution block"""
    print("Running mass balance gas volume verification plot...")
    error = plot_mass_balance('mass_balance_gas_volume', '')
    print("Mass balance gas volume verification completed successfully!")
