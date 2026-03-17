import numpy as np
import matplotlib.pyplot as plt
from matplotlib.transforms import blended_transform_factory
import fdsplotlib

plot_style = fdsplotlib.get_plot_style("fds")

nu = 1e-5
g  = 9.8

dx = np.logspace(-5, 2, 100)
plot_x_min = 1e-5
plot_x_max = 1e2
plot_y_min = 3e-6
plot_y_max = 3

k_sgs = 1.5 * (dx / np.pi)**(2/3)
tau_d = dx**2 / nu
tau_u = 0.75 * dx**(2/3)
tau_g = np.sqrt(2 * dx / g)
tau_chem = 1e-5*np.ones(len(dx)) # see Fox, p. 153
tau_flame = 1*np.ones(len(dx)) # puts a limit on the time scale when the flame volume is less than the cell volume
tau_min = np.maximum(dx / 100, tau_chem)

fig = fdsplotlib.plot_to_fig(x_data=dx,y_data=tau_d, marker_style='k--',linewidth=1,
                             plot_type='loglog', x_min=plot_x_min, x_max=plot_x_max, y_min=plot_y_min, y_max=plot_y_max,
                             x_label=r'$\Delta$ (m)', y_label=r'$\tau_{\rm mix}$ (s)')

fdsplotlib.plot_to_fig(x_data=dx, y_data=tau_u,   marker_style='k--', figure_handle=fig)
fdsplotlib.plot_to_fig(x_data=dx, y_data=tau_g,   marker_style='k--', figure_handle=fig)
fdsplotlib.plot_to_fig(x_data=dx, y_data=tau_chem,   marker_style='k--', figure_handle=fig)
fdsplotlib.plot_to_fig(x_data=dx, y_data=tau_flame,   marker_style='k--', figure_handle=fig)
fdsplotlib.plot_to_fig(x_data=dx, y_data=tau_min, marker_style='k--', linewidth=3, figure_handle=fig)

# plot actual model segments
dx1 = np.min(dx)
dx2 = np.sqrt(nu * tau_chem[0])
tau1 = tau_chem
tau2 = tau_chem
fdsplotlib.plot_to_fig(x_data=[dx1, dx2], y_data=[tau1, tau2], marker_style='k-', linewidth=3, figure_handle=fig)

dx1 = np.sqrt(nu * tau_chem[0])
dx2 = (0.75 * nu)**(0.75)
tau1 = tau_chem[0]
tau2 = dx2**2 / nu
dx_eta = dx2
tau_eta = tau2
fdsplotlib.plot_to_fig(x_data=[dx1, dx2], y_data=[tau1, tau2], marker_style='k-', linewidth=3, figure_handle=fig)

dx1 = (0.75 * nu)**(0.75)
dx2 = ((4/3) * np.sqrt(2 / g))**6
tau1 = dx1**2 / nu
tau2 = np.sqrt(2 * dx2 / g)
fdsplotlib.plot_to_fig(x_data=[dx1, dx2], y_data=[tau1, tau2], marker_style='k-', linewidth=3, figure_handle=fig)

dx1 = ((4/3) * np.sqrt(2 / g))**6
dx2 = 0.5*g*tau_flame[0]**2
tau1 = np.sqrt(2 * dx1 / g)
tau2 = np.sqrt(2 * dx2 / g)
fdsplotlib.plot_to_fig(x_data=[dx1, dx2], y_data=[tau1, tau2], marker_style='k-', linewidth=3, figure_handle=fig)

dx1 = 0.5*g*tau_flame[0]**2
dx2 = 1.5e2
tau1 = tau_flame[0]
tau2 = tau_flame[-1]
fdsplotlib.plot_to_fig(x_data=[dx1, dx2], y_data=[tau1, tau2], marker_style='k-', linewidth=3, figure_handle=fig)

ax = fig.axes[0]

ax.annotate(r'$\tau_{\rm min} = \Delta/100$',
            xy=(.63,.55), xycoords='figure fraction',
            xytext=(0.63, 0.43), textcoords='figure fraction',
            arrowprops=dict(arrowstyle="->", color='black'),
            fontname=plot_style['Font_Name'], fontsize=plot_style['Label_Font_Size'])

ax.annotate(r'$\tau_{\rm g} \sim \Delta^{1/2}$',
            xy=(0.65, 0.81), xycoords='figure fraction',
            xytext=(0.60, 0.68), textcoords='figure fraction',
            arrowprops=dict(arrowstyle="->", color='black'),
            fontname=plot_style['Font_Name'], fontsize=plot_style['Label_Font_Size'])

ax.scatter([dx_eta], [tau_eta],
           s=4000,
           facecolors='none',
           edgecolors='black',
           linewidths=1.0,
           linestyle='--')

ax.annotate(r'$\Delta = O(\eta)$',
            xy=(0.27, 0.61), xycoords='figure fraction',
            xytext=(0.20, 0.71), textcoords='figure fraction',
            arrowprops=dict(arrowstyle="->", color='black'),
            fontname=plot_style['Font_Name'], fontsize=plot_style['Label_Font_Size'])

ax.annotate(r'$\tau_{\rm u} \sim \Delta^{2/3}$',
            xy=(0.38, 0.60), xycoords='figure fraction',
            xytext=(0.39, 0.49), textcoords='figure fraction',
            arrowprops=dict(arrowstyle="->", color='black'),
            fontname=plot_style['Font_Name'], fontsize=plot_style['Label_Font_Size'])

ax.annotate(r'$\tau_{\rm d} \sim \Delta^2$',
            xy=(0.23, 0.38), xycoords='figure fraction',
            xytext=(0.27, 0.27), textcoords='figure fraction',
            arrowprops=dict(arrowstyle="->", color='black'),
            fontname=plot_style['Font_Name'], fontsize=plot_style['Label_Font_Size'])

ax = fig.axes[0]
ax.set_position([0.16, 0.14, 0.73, 0.8])

trans = blended_transform_factory(ax.transAxes, ax.transData)

ax.text(1.02, tau_chem[0], r'$\tau_{\rm chem}$',
        transform=trans, clip_on=False,
        ha='left', va='center',
        fontname=plot_style['Font_Name'],
        fontsize=plot_style['Label_Font_Size'])

ax.text(1.02, tau_flame[0], r'$\tau_{\rm flame}$',
        transform=trans, clip_on=False,
        ha='left', va='center',
        fontname=plot_style['Font_Name'],
        fontsize=plot_style['Label_Font_Size'])

fig.savefig('reaction_time_scale.pdf', format='pdf')


