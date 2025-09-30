
import numpy as np
import matplotlib.pyplot as plt
import os
import fdsplotlib

plot_style = fdsplotlib.get_plot_style('fds')

expdir = '../../../exp/BRE_Spray/'
outdir = '../../../out/BRE_Spray/'
pltdir = '../../Manuals/FDS_Validation_Guide/SCRIPT_FIGURES/BRE_LEMTA_Spray/'

git_file = outdir + 'BRE_Spray_A_1_git.txt'
version_string = fdsplotlib.get_version_string(git_file)

# Load experimental data and FDS prediction
exp_data = np.zeros((29, 9))  # Initialize with larger size to accommodate LEMTA data later
exp_data[0:24, 0:7] = np.loadtxt(expdir + 'BRE_Spray_Test.csv', delimiter=',', skiprows=2)

# Nozzle names
Nozzle = ['A', 'B', 'D']

# Experimental markers and colors
marker = ['bd','gs','ro','kv']
color  = ['b','g','r','k']

# Experimental data row indices
exp_rows = np.array([[0, 1, 2, 3, 4, 5, 6, 7],
                     [8, 9, 10, 11, 12, 13, 14, 15],
                     [16, 17, 18, 19, 20, 21, 22, 23]])

# Initialize arrays for FDS data
FDS_Flux_0 = np.zeros((4, 8))  # 4 to accommodate LEMTA data later
FDS_d32 = np.zeros((4, 8))
FDS_w = np.zeros((4, 8))
FDS_Flux = np.zeros((4, 8))
FDS_Attenuation = np.zeros((4, 8))

# Set up plots

figd = fdsplotlib.plot_to_fig(x_data=[0,0], y_data=[0,0],
                              x_min=0, x_max=9, y_min=0, y_max=800,
                              revision_label=version_string,
                              x_label='Pressure (bar)', y_label='Mean Diameter ($\mu$m)',
                              legend_location='upper left')

figw = fdsplotlib.plot_to_fig(x_data=[0,0], y_data=[0,0],
                              x_min=0, x_max=9, y_min=0, y_max=8,
                              revision_label=version_string,
                              x_label='Pressure (bar)', y_label='Mean W-Velocity (m/s)',
                              legend_location='upper left')

figa = fdsplotlib.plot_to_fig(x_data=[0,40], y_data=[0,40], marker_style='k-',
                              x_min=0, x_max=40, y_min=0, y_max=40,
                              figure_size=(plot_style['Scat_Paper_Width'],plot_style['Scat_Paper_Height']),
                              plot_size=(plot_style['Scat_Plot_Width'],plot_style['Scat_Plot_Height']),
                              plot_origin=(plot_style['Scat_Plot_X'],plot_style['Scat_Plot_Y']),
                              revision_label=version_string,
                              x_label='Measured Attenuation (\%)',
                              y_label='Predicted Attenuation (\%)')

# Plot diameter and velocity for Nozzles A, B, and D

for n in range(3):  # 0, 1, 2 for nozzles A, B, D
    for p in range(8):  # 0 to 7 for pressure points
        
        FDS_File = outdir + 'BRE_Spray_' + Nozzle[n] + '_' + str(p+1) + '_devc.csv'
        fds_data = np.loadtxt(FDS_File, delimiter=',', skiprows=2)
        
        # initial flux between t = 0.1 and 1.0 s
        i1 = np.where(fds_data[:, 0] >= 0.1)[0][0]
        i2 = np.where(fds_data[:, 0] <= 4.0)[0][-1]
        FDS_Flux_0[n, p] = np.mean(fds_data[i1:i2+1, 1])
        
        # spray properties at t = 9 s.
        i2 = np.where(fds_data[:, 0] < 9)[0][-1] - 1
        FDS_d32[n, p] = fds_data[i2, 2] * 1E6
        FDS_w[n, p] = -1 * fds_data[i2, 3]
        
        # Final flux as mean between t=11 and 15 s.
        i1 = np.where(fds_data[:, 0] >= 12)[0][0]
        i2 = np.where(fds_data[:, 0] > 10)[0][-1] - 1
        FDS_Flux[n, p] = np.mean(fds_data[i1:i2+1, 1])
        FDS_Attenuation[n, p] = 100 * (FDS_Flux_0[n, p] - FDS_Flux[n, p]) / FDS_Flux_0[n, p]
    
    # plot dv50 points

    fdsplotlib.plot_to_fig(x_data=range(1,9), y_data=exp_data[exp_rows[n, :], 5], data_label=f'Exp Nozzle {Nozzle[n]}', 
                           marker_style=marker[n], marker_fill_color='none', figure_handle=figd)
    fdsplotlib.plot_to_fig(x_data=range(1,9), y_data=FDS_d32[n, :], data_label=f'FDS Nozzle {Nozzle[n]}', 
                           marker_style=marker[n], marker_fill_color=color[n],figure_handle=figd)

    # plot w points

    fdsplotlib.plot_to_fig(x_data=range(1,9), y_data=exp_data[exp_rows[n, :], 6], data_label=f'Exp Nozzle {Nozzle[n]}', 
                           marker_style=marker[n], marker_fill_color='none', figure_handle=figw)
    fdsplotlib.plot_to_fig(x_data=range(1,9), y_data=FDS_w[n, :], data_label=f'FDS Nozzle {Nozzle[n]}', 
                           marker_style=marker[n], marker_fill_color=color[n], figure_handle=figw)

    # plot attenuation points for Nozzles A, B and D

    fdsplotlib.plot_to_fig(x_data=exp_data[exp_rows[n, :],3], y_data=FDS_Attenuation[n, :],
                           marker_style=marker[n], figure_handle=figa, data_label=f'BRE (Nozzle {Nozzle[n]})')

# Save velocity and diameter plots

figd.savefig(pltdir + 'BRE_Spray_Diameter.pdf', format='pdf')
plt.close(figd)
figw.savefig(pltdir + 'BRE_Spray_W.pdf', format='pdf')
plt.close(figw)

# Add LEMTA attenuation data to BRE data

expdir = '../../../exp/LEMTA_Spray/'
outdir = '../../../out/LEMTA_Spray/'

# load experimental data and FDS prediction
exp_data[0:5, 0:4] = np.loadtxt(expdir + 'LEMTA_Spray_Test.csv', delimiter=',', skiprows=2, usecols=(0, 1, 2, 3))

n = 3  # Index for LEMTA data
for p in range(5):
    
    FDS_File = outdir + 'LEMTA_Spray_' + str(p+1) + '_devc.csv'
    fds_data = np.loadtxt(FDS_File, delimiter=',', skiprows=2)
    
    # initial flux between t = 1 and 10 s
    i1 = np.where(fds_data[:, 0] > 1.0)[0][0]
    i2 = np.where(fds_data[:, 0] <= 10.0)[0][-1]
    FDS_Flux_0[n, p] = np.mean(fds_data[i1:i2+1, 1])
    
    # Final flux as mean between t=10 and 15 s
    i1 = np.where(fds_data[:, 0] > 10)[0][0]
    i2 = np.where(fds_data[:, 0] <= 30)[0][-1] - 1
    FDS_Flux[n, p] = np.mean(fds_data[i1:i2+1, 1])
    FDS_Attenuation[n, p] = 100 * (FDS_Flux_0[n, p] - FDS_Flux[n, p]) / FDS_Flux_0[n, p]

fdsplotlib.plot_to_fig(x_data=exp_data[0:5,3], y_data=FDS_Attenuation[n,0:5],
                       marker_style=marker[3], figure_handle=figa, data_label='LEMTA')

figa.savefig(pltdir + 'BRE_LEMTA_Spray_Attenuation.pdf', format='pdf')
plt.close(figa)


