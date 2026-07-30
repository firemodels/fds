
import pandas as pd
import numpy as np
import os

# include FDS plot styles, etc.
import fdsplotlib

# Get plot style parameters
plot_style = fdsplotlib.get_plot_style('fds')

expdir = '../../../exp/Sandia_Fireballs/'
outdir = '../../../out/Sandia_Fireballs/'
figdir = '../../Manuals/FDS_Validation_Guide/SCRIPT_FIGURES/Sandia_Fireballs/'

E_file = ['ethane_fireball.csv','ethylene_fireball.csv','isopentane_fireball.csv']
M_file = ['ethane_fireball_cat_devc.csv','ethylene_fireball_cat_devc.csv','isopentane_fireball_cat_devc.csv']

marker = ['ko','ro','go']
label  = ['Ethane','Ethylene','Isopentane']

git_file =  os.path.join(outdir, f'ethane_fireball_cat_git.txt')
version_string = fdsplotlib.get_version_string(git_file)

fig = fdsplotlib.plot_to_fig([0,80],[0,80], marker_style='k-',
                             x_min=0, x_max=80, y_min=0, y_max=80,
                             figure_size=(plot_style['Scat_Paper_Width'],plot_style['Scat_Paper_Height']),
                             plot_size=(plot_style['Scat_Plot_Width'],plot_style['Scat_Plot_Height']),
                             plot_origin=(plot_style['Scat_Plot_X'],plot_style['Scat_Plot_Y']),
                             legend_location='upper left',
                             revision_label=version_string,
                             x_label='Measured Heat Flux (kW/m²)',
                             y_label='Predicted Heat Flux (kW/m²)'
                             )

for i in range(3):

    # Read first CSV file and extract second column into array E
    df1 = pd.read_csv(expdir + E_file[i], skiprows=0)
    E = df1.iloc[:,1].values  # Second column (index 1)
    
    # Read second CSV file, skip first two rows, and extract the max value into array M
    df2 = pd.read_csv(outdir + M_file[i], skiprows=2)
    M = df2.iloc[:, 6:].max().values
    
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
                                 
    df = pd.DataFrame({'Index': range(len(M_clean)),'Heat Flux': M_clean})
    df.to_csv(outdir + label[i] + '.csv', index=False)
    df = pd.DataFrame({'Index': range(len(E_clean)),'Heat Flux': E_clean})
    df.to_csv(expdir + label[i] + '.csv', index=False)

# Save as PDF
fig_file =  os.path.join(figdir, f'Sandia_Fireballs.pdf')
fig.savefig(fig_file, format='pdf')

