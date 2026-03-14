import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
import fdsplotlib

plot_style = fdsplotlib.get_plot_style("fds")

nu = 1e-5
g  = 9.8

dx = np.logspace(-5, 2, 100)

k_sgs = 1.5 * (dx / np.pi)**(2/3)
tau_d = dx**2 / nu
tau_u = 0.75 * dx**(2/3)
tau_g = np.sqrt(2 * dx / g)
tau_chem = 1e-5 # see Fox, p. 153
tau_min = np.maximum(dx / 100, tau_chem)

fig = fdsplotlib.plot_to_fig(x_data=dx,y_data=tau_d, marker_style='k--',linewidth=1,
                             plot_type='loglog', x_min=1e-5, x_max=1e2, y_min=1e-5, y_max=10,
                             x_label=r'$\Delta$ (m)', y_label=r'$\tau_{\rm mix}$ (s)')

fdsplotlib.plot_to_fig(x_data=dx, y_data=tau_u,   marker_style='k--', figure_handle=fig)
fdsplotlib.plot_to_fig(x_data=dx, y_data=tau_g,   marker_style='k--', figure_handle=fig)
fdsplotlib.plot_to_fig(x_data=dx, y_data=tau_min, marker_style='k--', linewidth=3, figure_handle=fig)

# plot actual model segments
dx1 = np.min(dx)
dx2 = np.sqrt(nu * tau_chem)
tau1 = tau_chem
tau2 = tau_chem
fdsplotlib.plot_to_fig(x_data=[dx1, dx2], y_data=[tau1, tau2], marker_style='k-', linewidth=3, figure_handle=fig)

dx1 = np.sqrt(nu * tau_chem)
dx2 = (0.75 * nu)**(0.75)
tau1 = tau_chem
tau2 = dx2**2 / nu
fdsplotlib.plot_to_fig(x_data=[dx1, dx2], y_data=[tau1, tau2], marker_style='k-', linewidth=3, figure_handle=fig)

dx1 = (0.75 * nu)**(0.75)
dx2 = ((4/3) * np.sqrt(2 / g))**6
tau1 = dx1**2 / nu
tau2 = np.sqrt(2 * dx2 / g)
fdsplotlib.plot_to_fig(x_data=[dx1, dx2], y_data=[tau1, tau2], marker_style='k-', linewidth=3, figure_handle=fig)

dx1 = ((4/3) * np.sqrt(2 / g))**6
dx2 = 100
tau1 = np.sqrt(2 * dx1 / g)
tau2 = np.sqrt(2 * dx2 / g)
fdsplotlib.plot_to_fig(x_data=[dx1, dx2], y_data=[tau1, tau2], marker_style='k-', linewidth=3, figure_handle=fig)

ax = fig.axes[0]

ax.annotate(r'$\tau_{\rm min} \sim \Delta$',
            xy=(4/7,1/3), xycoords='figure fraction',
            xytext=(0.65, 0.25), textcoords='figure fraction',
            arrowprops=dict(arrowstyle="->", color='black'),
            fontname=plot_style['Font_Name'], fontsize=plot_style['Label_Font_Size'])

ax.annotate(r'$\tau_{\rm g} \sim \Delta^{1/2}$',
            xy=(0.65, 0.70), xycoords='figure fraction',
            xytext=(0.70, 0.65), textcoords='figure fraction',
            arrowprops=dict(arrowstyle="->", color='black'),
            fontname=plot_style['Font_Name'], fontsize=plot_style['Label_Font_Size'])

ellipse = Ellipse((0.23 + 0.15/2, 0.35 + 0.2/2), 0.15, 0.2,
                  edgecolor='black', facecolor='none',
                  linestyle='--', linewidth=1, transform=fig.transFigure)
ax.add_patch(ellipse)

ax.annotate(r'$\Delta = O(\eta)$',
            xy=(0.28, 0.55), xycoords='figure fraction',
            xytext=(0.25, 0.65), textcoords='figure fraction',
            arrowprops=dict(arrowstyle="->", color='black'),
            fontname=plot_style['Font_Name'], fontsize=plot_style['Label_Font_Size'])

ax.annotate(r'$\tau_{\rm u} \sim \Delta^{2/3}$',
            xy=(0.41, 0.535), xycoords='figure fraction',
            xytext=(0.445, 0.46), textcoords='figure fraction',
            arrowprops=dict(arrowstyle="->", color='black'),
            fontname=plot_style['Font_Name'], fontsize=plot_style['Label_Font_Size'])

ax.annotate(r'$\tau_{\rm d} \sim \Delta^2$',
            xy=(0.25, 0.32), xycoords='figure fraction',
            xytext=(0.30, 0.29), textcoords='figure fraction',
            arrowprops=dict(arrowstyle="->", color='black'),
            fontname=plot_style['Font_Name'], fontsize=plot_style['Label_Font_Size'])

fig.savefig('reaction_time_scale.pdf', format='pdf')


