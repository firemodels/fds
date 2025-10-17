
# Spray pattern plots are saved to fds/Manuals/FDS_User_Guide/FIGURES/spray_pattern_mu...pdf

import numpy as np
import matplotlib.pyplot as plt
import fdsplotlib

plot_style = fdsplotlib.get_plot_style('fds')

pltdir = '../../Manuals/FDS_User_Guide/FIGURES/'

b = [0, 5, 10, 100, 1000]
mu = np.array([0, 10, 20]) * np.pi / 180
mu_str = ['0', '10', '20']

phi_min = 0
phi_max = np.array([22.5, 45, 90]) * np.pi / 180
phi_str = ['22p5', '45', '90']
mark = ['k-','r-','g-','m-','c-']

for k in range(len(phi_max)):
    phi = np.linspace(phi_min, phi_max[k], 1000)
    for j in range(len(mu)):
        fig = fdsplotlib.plot_to_fig(x_data=[-1,-1], y_data=[-1,-1],
                             x_min=0, x_max=90, y_min=0, y_max=1.1,
                             x_label=r'$\phi$ (degrees)',
                             y_label=r'$f(\phi)$')

        for i in range(len(b)):
            f = np.exp(-b[i] * ((phi - mu[j]) / (phi_max[k] - phi_min)) ** 2)
            fdsplotlib.plot_to_fig(x_data=phi*180/np.pi, y_data=f, marker_style=mark[i], data_label=r'$\beta=$' + str(b[i]), figure_handle=fig)

        ax = plt.gca()
        ax.invert_yaxis()
        ax.text(45, 0.1, r'$\mu=$' + str(mu[j]*180/np.pi) + r', $\phi_{\rm max}=$' + str(phi_max[k]*180/np.pi) + ' degrees', fontsize=plot_style['Key_Font_Size'])

        plt.savefig(pltdir + 'spray_pattern_mu_' + mu_str[j] + '_phimax_' + phi_str[k] + '.pdf', format='pdf')
        plt.close()

