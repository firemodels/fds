
# This script creates additional figures related to the USFS_Catchpole cases.
# The objective is to better identify trends/issues within the parameter space.

import numpy as np
import matplotlib.pyplot as plt
import os, sys
import pandas as pd
import matplotlib as mpl

filedir = os.path.dirname(__file__)
firemodels = os.path.join(filedir,'..','..','..','..')
sys.path.append(filedir+os.sep+'..'+os.sep)
import fdsplotlib

plot_style = fdsplotlib.get_plot_style("fds")

base_path = os.path.join(firemodels,'out','USFS_Catchpole')
fig_path = os.path.join(firemodels,'fds','Manuals','FDS_Validation_Guide','SCRIPT_FIGURES','USFS_Catchpole')
validation_path = os.path.join(firemodels,'fds','Validation','USFS_Catchpole','FDS_Input_Files')

# Experiment parameters
tests = pd.read_csv(os.path.join(validation_path,"Test_Matrix.csv"))

# Calculate spread rates for summary plotting
for ti,test in tests.iterrows():
    R = tests['R'].iloc[ti]
    chid = tests['Test'].iloc[ti]
    fds_file = os.path.join(base_path, f"{chid}_devc.csv")
    git_file = os.path.join(base_path, f"{chid}_git.txt")
    version_string = fdsplotlib.get_version_string(git_file)
    
    fds_data = pd.read_csv(fds_file,header=1)
    
    # fit spread rate
    try:
        # fit slope filtering to positions between 2 m and 7 m from ignition
        xi = (fds_data['x']>2) & (fds_data['x']<7)
        R_FDS,intercept = np.polyfit(fds_data[xi]['Time'], fds_data[xi]['x'], 1)
    # not enough data to fit 
    except:
        R_FDS=0.
        
    if R_FDS<0:
        R_FDS=0
    
    # add fds data to full table for summary plotting      
    tests.loc[ti,'R_FDS'] = R_FDS
    # No spread if zero from above or the fire does not reach at least 7 m
    if (R_FDS < 1e-5) or (fds_data['x'].max() < 7):
        tests.loc[ti,'category'] = 'no spread'
    else:
        tests.loc[ti,'category'] = 'spread'

    
# Create summary plots

dep_variables={"s":"Surface-to-Volume Ratio, s (1/m)",
               "beta":"Packing Ratio, beta (-)",
               "U":"Wind Speed, U (m/s)",
               "M":"Moisture Content, M (-)"}

fuel_labels=["MF","PPMC","EXSC","EX"]
fuel_names=["Pine sticks","Pine needles","Coarse excelsior","Regular excelsior"]
colors = ['b','g','r','c','m','y','k'] # matlab defauls

for dvar in dep_variables:
    fig_file = os.path.join(fig_path, f"Catchpole_R_v_{dvar}.pdf")
    
    # show +/- 20% relative error
    [xmin,xmax] = [0.8*tests[dvar].min(),1.1*tests[dvar].max()]
    fig = fdsplotlib.plot_to_fig(x_data=[xmin, xmax], y_data=[0.8,0.8],
                                 plot_type='semilogy', marker_style='k--',
                                 x_min=xmin, x_max=xmax, y_min=0.1, y_max=15,
                                 x_label=dep_variables[dvar], y_label=r"$\mathrm{R_{FDS}/R_{Exp}}$",
                                 revision_label=version_string)

    fdsplotlib.plot_to_fig(x_data=[xmin, xmax], y_data=[1.2,1.2], marker_style='k--', figure_handle=fig)

    for i in range(0, len(fuel_labels)):
        fuel = fuel_labels[i]
        filtered_data = tests[tests['Test'].str.startswith(fuel)]
        if fuel=='EX':
            filtered_data = tests[
                (tests['Test'].str.startswith(fuel))&(~tests['Test'].str.startswith('EXSC'))]
        fdsplotlib.plot_to_fig(x_data=filtered_data[dvar], y_data=filtered_data['R_FDS']/filtered_data['R'],
                               data_label=fuel_names[i], marker_style=colors[i]+'o', figure_handle=fig)

    plt.savefig(fig_file)
    plt.close(fig)

# plot no-spread conditions

go_mask    = tests['category'] == 'spread'
nogo_mask  = tests['category'] == 'no spread'

# Create scatter plot for each fuel type
for i, fuel in enumerate(fuel_labels):
    fuel_name = fuel_names[i]
    fig_file = os.path.join(fig_path, f"Catchpole_no_spread_{fuel}.pdf")
    
    # Filter for this fuel type
    if fuel == 'EX':
        fuel_mask = tests['Test'].str.startswith(fuel) & (~tests['Test'].str.startswith('EXSC'))
    else:
        fuel_mask = tests['Test'].str.startswith(fuel)
    
    # Get data for this fuel type
    fuel_go_data = tests[go_mask & fuel_mask]
    fuel_nogo_data = tests[nogo_mask & fuel_mask]
    
    s_value = fuel_go_data['s'].iloc[0]
    # Create figure
    fig = fdsplotlib.plot_to_fig(x_data=[-1,1], y_data=[0,0],
                                    marker_style='k--',
                                    x_label='Moisture Content, M (-)',
                                    y_label='Wind Speed, U (m/s)',
                                    x_min=0,
                                    x_max=0.3,
                                    y_min=-0.2,
                                    y_max=3.5,
                                    revision_label=version_string,
                                    plot_title=rf'{fuel_name}, {s_value:.0f} m$^{{-1}}$')
    ax = plt.gca()
    
    # Fixed scale for beta (0 to 0.1)
    beta_min_fixed = 0.0
    beta_max_fixed = 0.1
    
    # Plot spread cases (background, gray, sized by beta)
    if len(fuel_go_data) > 0:
        sizes_go = (fuel_go_data['beta'] - beta_min_fixed) / (beta_max_fixed - beta_min_fixed) * 100 + 20
        ax.scatter(fuel_go_data['M'], fuel_go_data['U'], 
                  s=sizes_go, c='gray', alpha=0.3, edgecolors='none')
    
    # Plot no-spread cases (foreground, colored by beta, sized by beta)
    if len(fuel_nogo_data) > 0:
        cmap = plt.cm.plasma
        norm = mpl.colors.Normalize(vmin=beta_min_fixed, vmax=beta_max_fixed)
        sizes_nogo = (fuel_nogo_data['beta'] - beta_min_fixed) / (beta_max_fixed - beta_min_fixed) * 100 + 20
        scatter = ax.scatter(fuel_nogo_data['M'], fuel_nogo_data['U'], 
                  s=sizes_nogo, c=fuel_nogo_data['beta'], cmap=cmap, norm=norm,
                  alpha=0.8, edgecolors='black', linewidths=0.5)
        
        # Add colorbar
        cbar = plt.colorbar(scatter, ax=ax)
        cbar.set_label('Packing Ratio, beta (-)', rotation=90, labelpad=15, fontsize=plot_style['Label_Font_Size'])
        cbar.ax.tick_params(labelsize=plot_style['Label_Font_Size'])

    plt.savefig(fig_file)
    plt.close()
