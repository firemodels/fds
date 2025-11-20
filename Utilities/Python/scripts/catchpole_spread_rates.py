
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

fig_file = os.path.join(fig_path, "Catchpole_no_spread.pdf")
# Dummy call to establish figure
fig = fdsplotlib.plot_to_fig(x_data=[-1,-1], y_data=[-1,-1],
                                x_label='Parameter',y_label='Normalized value',
                                x_min=0,x_max=3,y_min=0,y_max=1,ynumticks=2,
                                revision_label=version_string)
ax = plt.gca()

# Normalize by max and min
tests_normalized = tests
tests_normalized[list(dep_variables.keys())] = tests[list(dep_variables.keys())].apply(
    lambda x: (x - x.min()) / (x.max() - x.min()))

go_mask    = tests_normalized['category'] == 'spread'
nogo_mask  = tests_normalized['category'] == 'no spread'

pd.plotting.parallel_coordinates(tests_normalized[go_mask], 'category', 
                                 cols=['s','beta','M','U'],
                                 color=[(0.,0.,0.,.04)],
                                 ax=ax,
                                 ls='-')
ax.get_legend().remove()  # Hide the parallel_coordinates legend

color_var = 'beta'   # or 'beta', 'M', 'U'
cmap = plt.cm.plasma 
norm = mpl.colors.Normalize(vmin=tests[color_var].min(),
                     vmax=tests[color_var].max())
colors = cmap(norm(tests[color_var]))

for idx, (i, row) in enumerate(tests_normalized[nogo_mask].iterrows()):
    x_vals = ['s','beta','M','U'] 
    y_vals = [row['s'], row['beta'], row['M'], row['U']]
    ax.plot(x_vals, y_vals, color=colors[i], linestyle='-')

plt.savefig(fig_file)
plt.close()
