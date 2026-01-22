
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import numpy as np
import pandas as pd
from scipy.interpolate import griddata
from fdsplotlib import configure_fds_fonts

configure_fds_fonts(usetex=True)

expdir = '../../../exp/Submodules/macfp-db/Wall_Fires/JIS_Facade/Experimental_Data/'
pltdir = '../../Manuals/FDS_Validation_Guide/SCRIPT_FIGURES/TUS_Facade/'

# Read the CSV file, skipping the header row
df = pd.read_csv(expdir + 'Sun_FAM_2024_mean_temperature.csv', skiprows=1, header=None)

# Extract coordinates
x = df.iloc[:, 0].values
z = df.iloc[:, 1].values

xi = np.linspace(np.min(x), np.max(x), 200)
zi = np.linspace(np.min(z), np.max(z), 200)
Xi, Zi = np.meshgrid(xi, zi)

custom_levels = np.linspace(48.0, 1175.0, 9)

# Make 3 sets of plots for 600, 750, 900 kW
for j in range(3):

    # Temperatures at the center, center-right and right
    T1 = df.iloc[:, 2+3*j].values
    T2 = df.iloc[:, 3+3*j].values
    T3 = df.iloc[:, 4+3*j].values
    
    # Create figure with three subplots
    fig, axes = plt.subplots(1, 3, figsize=(4.5, 4.5))
    fig.suptitle(f'{600+j*150} kW', fontsize=14)
    fig.subplots_adjust(left=0.02, right=0.85, wspace=-0.30)
    
    # Plot contour maps
    contours = []
    for i, (ax, T) in enumerate(zip(axes, [T1, T2, T3])):
        Ti = griddata((x, z), T, (Xi, Zi), method='cubic')
        if Ti.ndim == 3:
            Ti = Ti[..., 0]
        cs = ax.contourf(Xi, Zi, Ti, levels=custom_levels, cmap='rainbow', vmin=50, vmax=1100)
        contours.append(cs)
        
        # Set equal aspect ratio for same scaling
        ax.set_aspect('equal')
        
        # Only add labels to left and bottom
        if i == 0:  # leftmost plot gets y-axis labels
            ax.set_ylabel('Height (m)')
            ax.set_title('Middle', fontsize=10)
            ax.set_xlabel('Distance (m)')
            ax.xaxis.set_major_locator(MultipleLocator(0.4))
        elif i == 1:
            ax.set_title('Middle-Right', fontsize=10)
            ax.set_xticks([])
            ax.set_yticks([])
        else:
            ax.set_title('Right', fontsize=10)
            ax.set_xticks([])
            ax.set_yticks([])
        
    # Adjust layout and add colorbar to the right
    #plt.tight_layout(pad=0.0, w_pad=0.0, h_pad=0.0)
    cbar = plt.colorbar(contours[0], ax=axes, orientation='vertical', 
                       fraction=0.05, pad=0.04)
    cbar.set_label('Temperature (°C)')
    
    plt.savefig(pltdir + 'JIS_facade_exp_contours_' + str(600+j*150) + '.pdf')


