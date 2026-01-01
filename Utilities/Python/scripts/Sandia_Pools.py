
import pandas as pd
import numpy as np
import os

# include FDS plot styles, etc.
import fdsplotlib

# Get plot style parameters
plot_style = fdsplotlib.get_plot_style('fds')

expdir = '../../../exp/Sandia_Jets_Pools_Fireballs/'
outdir = '../../../out/Sandia_Jets_Pools_Fireballs/'
figdir = '../../Manuals/FDS_Validation_Guide/SCRIPT_FIGURES/Sandia_Jets_Pools_Fireballs/'

E_file = ['ethane_pool.csv','ethylene_pool.csv','propane_pool.csv','isopentane_pool.csv']
M_file = ['ethane_pool_cat_devc.csv','ethylene_pool_cat_devc.csv','propane_pool_cat_devc.csv','isopentane_pool_cat_devc.csv']

marker = ['ko','ro','go','bo']
label  = ['Ethane','Ethylene','Propane','Isopentane']

git_file =  os.path.join(outdir, f'ethane_pool_cat_git.txt')
version_string = fdsplotlib.get_version_string(git_file)

fig = fdsplotlib.plot_to_fig([0,50],[0,50], marker_style='k-',
                             x_min=0, x_max=50, y_min=0, y_max=50,
                             figure_size=(plot_style['Scat_Paper_Width'],plot_style['Scat_Paper_Height']),
                             plot_size=(plot_style['Scat_Plot_Width'],plot_style['Scat_Plot_Height']),
                             plot_origin=(plot_style['Scat_Plot_X'],plot_style['Scat_Plot_Y']),
                             legend_location='upper left',
                             revision_label=version_string,
                             x_label='Measured Heat Flux (kW/m²)',
                             y_label='Predicted Heat Flux (kW/m²)'
                             )

for i in range(4):

    # Read first CSV file and extract second column into array E
    df1 = pd.read_csv(expdir + E_file[i], skiprows=0)
    E = df1.iloc[:,1].values  # Second column (index 1)
    #col1 = df1.iloc[:,1].values  # Second column (index 1)
    #col2 = df1.iloc[:,2].values  # Third column (index 2)
    #E = np.concatenate([col1, col2])  # Combine both columns into single array
    
    # Read second CSV file, skip first two rows, and extract last row into array M
    df2 = pd.read_csv(outdir + M_file[i], skiprows=2)
    M = df2.iloc[-1,1:].values  # Last row, all columns
    
    # Remove NaN values by creating a mask for valid (non-NaN) values in both arrays
    # Ensure both arrays have the same length by taking the minimum length
    min_length = min(len(E), len(M))
    E = E[:min_length]
    M = M[:min_length]
    
    # Create mask for non-NaN values in both arrays
    valid_mask = ~(np.isnan(E) | np.isnan(M))
    E_clean = E[valid_mask]
    M_clean = M[valid_mask]
    
    # Add data to scatterplot
    fdsplotlib.plot_to_fig(E_clean, M_clean,figure_handle=fig,
                           marker_style=marker[i],
                           data_label=label[i])
                                 
# Save as PDF
fig_file =  os.path.join(figdir, f'Sandia_Pools.pdf')
fig.savefig(fig_file, format='pdf')

