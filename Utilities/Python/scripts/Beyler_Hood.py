
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import fdsplotlib

plot_style = fdsplotlib.get_plot_style('fds')

expdir = '../../../exp/Beyler_Hood/'
outdir = '../../../out/Beyler_Hood/'
pltdir = '../../Manuals/FDS_Validation_Guide/SCRIPT_FIGURES/Beyler_Hood/'

git_file = outdir + 'Beyler_Hood_acetone_117_lr_git.txt'
version_string = fdsplotlib.get_version_string(git_file)

# load experimental data and FDS prediction
exp_data = pd.read_csv(expdir + 'Beyler_Hood_data_lr.csv', skiprows=2, header=None).values

N_Fuels = 7
N_Species = 6

Fuel = ['acetone','ethanol','isopropanol','methanol','propane','propylene','toluene']
marker = ['ko','rs','bd','g^','mv','c>','k<'] # Marker for each fuel
NumPoints = [5,5,5,5,15,7,5]  # Number of points for each fuel

Species = ['O$_2$','CO$_2$','H$_2$O','CO','UHC','Soot']
SaveName = ['O2','CO2','H2O','CO','UHC','Soot']
Xmax = [0.25, 0.25, 0.2, 0.10, 0.10, 0.02]

TestID = np.zeros((N_Fuels,max(NumPoints)), dtype=int)  # Make it large enough to accommodate all data
TestID[0,0:NumPoints[0]] = [117,119,122,142,145]
TestID[1,0:NumPoints[1]] = [106,107,108,110,115]
TestID[2,0:NumPoints[2]] = [130,132,133,136,141]
TestID[3,0:NumPoints[3]] = [942,943,945,947,951]
TestID[4,0:NumPoints[4]] = [232,257,287,303,307,318,322,334,355,359,371,389,429,433,445]
TestID[5,0:NumPoints[5]] = [780,805,859,870,882,886,910]
TestID[6,0:NumPoints[6]] = [160,162,165,166,170]

ExpPlot = np.zeros((N_Fuels,max(NumPoints),N_Species))
FDSPlot = np.zeros((N_Fuels,max(NumPoints),N_Species))

for f in range(N_Fuels):
    for s in range(NumPoints[f]):
        FDS_File = outdir + 'Beyler_Hood_' + Fuel[f] + '_' + str(TestID[f, s]) + '_lr_devc.csv'
        fds_data = pd.read_csv(FDS_File, skiprows=2, header=None).values
        n_fds = fds_data.shape[0]
        for ns in range(N_Species):
            ExpPlot[f, s, ns] = exp_data[s, f * N_Species + ns]
            FDSPlot[f, s, ns] = np.mean(fds_data[n_fds-60:n_fds,ns+1])

for ns in range(N_Species):
    fig = fdsplotlib.plot_to_fig(x_data=[0,Xmax[ns]], y_data=[0,Xmax[ns]], marker_style='k-',
                                 x_min=0, x_max=Xmax[ns], y_min=0, y_max=Xmax[ns],
                                 figure_size=(plot_style['Scat_Paper_Width'],plot_style['Scat_Paper_Height']),
                                 plot_size=(plot_style['Scat_Plot_Width'],plot_style['Scat_Plot_Height']),
                                 plot_origin=(plot_style['Scat_Plot_X'],plot_style['Scat_Plot_Y']),
                                 revision_label=version_string,
                                 x_label='Measured ' + Species[ns] + ' (Mass Fraction)',
                                 y_label='Predicted ' + Species[ns] + ' (Mass Fraction)')

    for f in range(N_Fuels):
        fdsplotlib.plot_to_fig(x_data=ExpPlot[f,:,ns], y_data=FDSPlot[f,:,ns],
                               marker_style=marker[f], figure_handle=fig, data_label=Fuel[f])

    plotname = pltdir + 'Beyler_Hood_' + SaveName[ns] + '.pdf'
    plt.savefig(plotname, format='pdf')
    plt.close()
    
