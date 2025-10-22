"""
mass_balance_reac.py
Validates mass conservation for species with chemical reactions.
Original MATLAB script by McDermott (10 July 2020)
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import fdsplotlib
import os

def plot_mass_balance(chid, title_text, mass_id, devc_id, y_min=None, y_max=None, ynumticks=8):
    
    plot_style = fdsplotlib.get_plot_style('fds')
    fds_dir = os.path.normpath(os.path.join(os.path.dirname(__file__),'..','..','..'))
    ddir = os.path.join(fds_dir, 'Verification','Species','')
    pltdir = os.path.join(fds_dir, 'Manuals','FDS_Verification_Guide','SCRIPT_FIGURES','')
    
    mass_file = os.path.join(ddir, f'{chid}_mass.csv')
    devc_file = os.path.join(ddir, f'{chid}_devc.csv')

    try:
        M = pd.read_csv(mass_file, skiprows=1)  # Skip units row, use column names as header
        F = pd.read_csv(devc_file, skiprows=1)  # Skip units row, use column names as header
    except FileNotFoundError as e:
        print(f'Error: File does not exist. Skipping case. {e}')
        return 0

    t = M.iloc[:, 0].values
    m = M[mass_id].values
    dmdt = np.zeros(len(t))
    for i in range(1, len(t)):
        dmdt[i] = (m[i] - m[i-1]) / (t[i] - t[i-1])
    
    mdot_out = (F[f'{devc_id} xmax'].values + 
                F[f'{devc_id} xmin'].values + 
                F[f'{devc_id} ymin'].values + 
                F[f'{devc_id} ymax'].values + 
                F[f'{devc_id} Burner'].values)
    
    # Generation/consumption from reaction
    gen = F[f'{devc_id} mdot reac'].values   
    # Balance: -dm/dt + (in - out) + generation (should be ~0)
    bal = -dmdt + mdot_out + gen
        
    # Get version string if git file exists
    git_file = os.path.join(ddir, f'{chid}_git.txt')
    version_string = fdsplotlib.get_version_string(git_file) if os.path.exists(git_file) else ''
    
    # Create plot with all data
    fig = fdsplotlib.plot_to_fig(x_data=t, y_data=dmdt, marker_style='g-',
                                 revision_label=version_string, legend_location='lower left', x_label='Time (s)', y_label='Mass Flow (kg/s)',
                                 data_label='accumulation',x_min=0, x_max=20, y_min=y_min, y_max=y_max, xnumticks=5, ynumticks=ynumticks)
    fdsplotlib.plot_to_fig(x_data=t, y_data=mdot_out, figure_handle=fig, marker_style='b-', data_label='in - out')
    fdsplotlib.plot_to_fig(x_data=t, y_data=gen, figure_handle=fig, marker_style='r-', data_label='generation')
    fdsplotlib.plot_to_fig(x_data=t, y_data=bal, figure_handle=fig, marker_style='k-', data_label='balance')
    ax = plt.gca()
    xl = ax.get_xlim()
    yl = ax.get_ylim()
    x_pos = xl[0] + 0.05 * (xl[1] - xl[0])
    y_pos = yl[0] + 0.9 * (yl[1] - yl[0])
    ax.text(x_pos, y_pos, title_text, fontsize=plot_style['Title_Font_Size'], 
            fontname=plot_style['Font_Name'])
    
    # Check balance and report error (normalized by total mass)
    mass_error = np.abs(np.mean(bal)) / m[-1] if m[-1] != 0 else np.nan
    if mass_error > 1.5e-3:
        print(f'Python Warning: mass error = {mass_error:.6e} in {chid} for species {mass_id}')
    
    # Save figure
    output_file = os.path.join(pltdir, f'{chid}_{devc_id}.pdf')
    plt.savefig(output_file, format='pdf')
    plt.close(fig)
        
    return mass_error

if __name__ == '__main__':
    """Main execution block"""    
    print("Running mass balance reac verification plots...")
    error1 = plot_mass_balance('mass_balance_reac', 'Propane Mass Balance', 'PROPANE', 'C3H8',y_min=-0.15,y_max=0.1,ynumticks=6)
    error2 = plot_mass_balance('mass_balance_reac', 'Oxygen Mass Balance', 'OXYGEN', 'O2',y_min=-2,y_max=1.5,ynumticks=8)
    error3 = plot_mass_balance('mass_balance_reac', 'Nitrogen Mass Balance', 'NITROGEN', 'N2',y_min=-6,y_max=4,ynumticks=6)
    error4 = plot_mass_balance('mass_balance_reac', 'Carbon Dioxide Mass Balance', 'CARBON DIOXIDE', 'CO2',y_min=-0.4,y_max=0.4,ynumticks=5)
    error5 = plot_mass_balance('mass_balance_reac', 'Water Vapor Mass Balance', 'WATER VAPOR', 'H2O',y_min=-0.2,y_max=0.2,ynumticks=5)
    error6 = plot_mass_balance('mass_balance_reac', 'Soot Mass Balance', 'SOOT', 'Soot',y_min=-0.002,y_max=0.002,ynumticks=5)
    print("Mass balance reac verification completed successfully!")
