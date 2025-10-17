
# Generate the Rosin-Rammler / log-normal particle size distribution plot in the FDS Tech Guide 

import numpy as np
import matplotlib.pyplot as plt
import fdsplotlib 

plot_style = fdsplotlib.get_plot_style("fds")

pltdir = '../../Manuals/FDS_Technical_Reference_Guide/SCRIPT_FIGURES/'

D_v50 = 1
gamma = 2.4
sigma = 1.15/gamma
N = 100
D = np.linspace(0, 3, N)
dD = D[1] - D[0]

F1 = np.zeros(N)
F2 = np.zeros(N)
Fv = np.zeros(N)
fv = np.zeros(N)
Fn = np.zeros(N)

for i in range(N):

    if i == 0:
        F1[i] = 0
    else:
        DP = 0.5 * (D[i-1] + D[i])
        F1[i] = F1[i-1] + np.exp(-np.log(DP/D_v50)**2/(2*sigma**2)) * dD / (np.sqrt(2*np.pi)*sigma*DP)

    F2[i] = 1 - np.exp(-0.693*(D[i]/D_v50)**gamma)

    if D[i] < D_v50:
        Fv[i] = F1[i]
    else:
        Fv[i] = F2[i]

    # probability density function (derivative of CDF)
    if i == 0:
        fv[i] = 0
    else:
        fv[i] = (Fv[i] - Fv[i-1]) / (D[i] - D[i-1])

    if i == 0:
        Fn[i] = 0
    else:
        DP = 0.5 * (D[i-1] + D[i])
        Fn[i] = Fn[i-1] + fv[i]/DP**3 * dD

Fn = Fn / Fn[-1]

fig = fdsplotlib.plot_to_fig(x_data=[D_v50, D_v50], y_data=[0, 1], marker_style='k:', data_label=r'$D_{\rm v,0.5}$',
                             x_min=0, x_max=3, y_min=0, y_max=1,
                             x_label='Droplet Diameter (mm)')

fdsplotlib.plot_to_fig(x_data=D, y_data=F1, marker_style='r-.', data_label=r'$F_{\rm v}$ Log-Normal', figure_handle=fig)
fdsplotlib.plot_to_fig(x_data=D, y_data=F2, marker_style='c-.', data_label=r'$F_{\rm v}$ Rosin-Rammler', figure_handle=fig)
fdsplotlib.plot_to_fig(x_data=D, y_data=Fv, marker_style='k-',  data_label=r'$F_{\rm v}$ Combined (FDS)', figure_handle=fig)
fdsplotlib.plot_to_fig(x_data=D, y_data=Fn, marker_style='b-',  data_label=r'$F_{\rm n}$ Cumulative Number Fraction', figure_handle=fig)

plt.savefig(pltdir + 'particle_size_distribution.pdf', format='pdf')
plt.close()

