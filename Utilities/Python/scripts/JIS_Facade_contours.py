
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import numpy as np
import pandas as pd
from scipy.interpolate import griddata

expdir = '../../../exp/Submodules/macfp-db/Wall_Fires/JIS_Facade/Experimental_Data/'
outdir = '../../../out/TUS_Facade/'
pltdir = '../../Manuals/FDS_Validation_Guide/SCRIPT_FIGURES/TUS_Facade/'

# Read the CSV file, skipping the first two header rows
df = pd.read_csv(outdir + 'JIS_facade_2cm_line.csv', skiprows=2, header=None)

# Remove rows with all NaN values and handle NaN values
df = df.dropna(how='all')

# Set vector z to column 5 (assuming 1-indexed, so column index 4)
z = df.iloc[:, 4].dropna().values

# Create vector x to be the distances from the wall
x = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8]

X, Z = np.meshgrid(x, z)

custom_levels = np.linspace(48.0, 1175.0, 9)

# Make 3 sets of plots for 600, 750, 900 kW
for j in range(3):

    # Temperatures at the center, center-right and right
    T1 = df.iloc[:, 5+27*j:14+27*j].dropna().values
    T2 = df.iloc[:,14+27*j:23+27*j].dropna().values
    T3 = df.iloc[:,23+27*j:32+27*j].dropna().values
    
    # Create figure with three subplots
    fig, axes = plt.subplots(1, 3, figsize=(6, 4.5))
    fig.suptitle(f'Temperature Contours, {600+j*150} kW', fontsize=14)
    fig.subplots_adjust(left=-0.05, right=0.90, wspace=-0.20)
    
    # Plot contour maps
    contours = []
    for i, (ax, T) in enumerate(zip(axes, [T1, T2, T3])):
        cs = ax.contourf(X, Z, T, levels=custom_levels, cmap='rainbow', vmin=50, vmax=1100)
        contours.append(cs)
        
        # Set equal aspect ratio for same scaling
        ax.set_aspect('equal')
        
        # Only add labels to left and bottom
        if i == 0:  # leftmost plot gets y-axis labels
            ax.set_ylabel('Height (m)')
            ax.set_title('Middle', fontsize=12)
        elif i == 1:
            ax.tick_params(labelleft=False)
            ax.set_title('Middle-Right', fontsize=12)
        else:
            ax.tick_params(labelleft=False)
            ax.set_title('Right', fontsize=12)
        
        # All plots get x-axis labels (bottom)
        ax.set_xlabel('Distance (m)')
        ax.xaxis.set_major_locator(MultipleLocator(0.4))
    
    # Adjust layout and add colorbar to the right
    #plt.tight_layout(pad=0.0, w_pad=0.0, h_pad=0.0)
    cbar = plt.colorbar(contours[0], ax=axes, orientation='vertical', 
                       fraction=0.05, pad=0.04)
    cbar.set_label('Temperature (°C)')
    
    plt.savefig(pltdir + 'JIS_facade_contours_' + str(600+j*150) + '.pdf')



# Read the experimental CSV file
#data = pd.read_csv(expdir + 'Sun_FAM_2024_mean_temperature.csv', skiprows=1)

# Extract x, z, T from first 3 columns
#x_exp  = data.iloc[:, 0].values
#z_exp  = data.iloc[:, 1].values
#T1_exp = data.iloc[:, 2].values
#T2_exp = data.iloc[:, 3].values
#T3_exp = data.iloc[:, 4].values

# Create regular grid for interpolation
#xi = np.linspace(x_exp.min(), x_exp.max(), 100)
#zi = np.linspace(z_exp.min(), z_exp.max(), 100)
#Xi, Zi = np.meshgrid(xi, zi)

# Interpolate unstructured data onto regular grid
#T1i_exp = griddata((x_exp, z_exp), T1_exp, (Xi, Zi), method='linear')
#T2i_exp = griddata((x_exp, z_exp), T2_exp, (Xi, Zi), method='linear')
#T3i_exp = griddata((x_exp, z_exp), T3_exp, (Xi, Zi), method='linear')

