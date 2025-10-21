# Toms
# 8-8-14
# backward_facing_step.m
# Converted by Floyd
#
# 10/20/2025

import numpy as np
import pandas as pd
import os
import math
import matplotlib.pyplot as plt
import fdsplotlib


def tight_subplot(Nh, Nw, gap=(0.01, 0.01), marg_h=(0.13, 0.08), marg_w=(0.11, 0.019)):
   '''
   Creates a grid of subplot axes with adjustable gaps and margins.
   Returns the Figure object and a flat list of Axes handles.
   '''

   if np.isscalar(gap):
      gap = (gap, gap)
   if np.isscalar(marg_w):
      marg_w = (marg_w, marg_w)
   if np.isscalar(marg_h):
      marg_h = (marg_h, marg_h)

   axh = (1.0 - sum(marg_h) - (Nh - 1) * gap[0]) / Nh
   axw = (1.0 - sum(marg_w) - (Nw - 1) * gap[1]) / Nw

   py = 1.0 - marg_h[1] - axh

   ha = []
   fig = plt.figure(figsize=(12, 8))

   for ih in range(Nh):
      px = marg_w[0]
      for ix in range(Nw):
         # Position: [left, bottom, width, height]
         pos = [px, py, axw, axh]
         ax = fig.add_axes(pos)
         ha.append(ax)
         px = px + axw + gap[1]
      py = py - axh - gap[0]

   return fig, ha

expdir = '../../../exp/Backward_Facing_Step/'
datdir = '../../../out/Backward_Facing_Step/'
pltdir = '../../Manuals/FDS_Validation_Guide/SCRIPT_FIGURES/Backward_Facing_Step/'
gitname= 'backward_facing_step_20_git.txt'
git_file = datdir+gitname
version_string = fdsplotlib.get_version_string(git_file)


rkappa = 1/.41
B = 5.2
rho = 1.1992662 # density, kg/m^3
U_0 = 7.72      # reference velocity, m/s
h = 0.0098      # step height, m
visc = 1.7801E-5
visc_nu = visc/rho
x_loc = ['-3','4','6','10']

nx = [5, 10, 20]
lnx = len(nx)
fds_marker = ['rs-','go-','bd-']
subcol = ['r','g','b']
subsym = ['s','o','d']
fds_key = ['$h/\\delta x$=5','$h/\\delta x$=10','$h/\\delta x$=20']

# --- File Existence Check ---
exp_file = os.path.join(expdir, 'backward_facing_step_data.csv')
if not os.path.exists(exp_file):
   print(f'Error: File {exp_file} does not exist. Skipping case.')
   quit()

data_files = []
for n in nx:
   dat_file = os.path.join(datdir, f'backward_facing_step_{n}_line.csv')
   if not os.path.exists(dat_file):
      print(f'Error: File {dat_file} does not exist. Skipping case.')
      quit()
   data_files.append(dat_file)

# --- Read Experimental and FDS Data ---

# Experimental Data (D): Assumed header is the first row (header=0)
D = pd.read_csv(exp_file)

# FDS Data (M): Assumed header is the second row (header=1)
M = []
for df_path in data_files:
   M.append(pd.read_csv(df_path, skiprows=1))


# --- 1. Streamwise Cf data along bottom of channel ---

# Exp Data (Cf-x/h)
xoh_cf = D['Cf-x/h'].values
Cf_exp = D['Cf'].values

error = 0.0005 * np.ones_like(Cf_exp)

fig = fdsplotlib.plot_to_fig(x_data=xoh_cf, y_data=Cf_exp, y_error=error,marker_style='ko',
      revision_label=version_string,x_min=0,x_max=20,y_min=-0.004,y_max=0.004,
      data_label='J&D',
      x_label=r'$x/h$',
      y_label=r'$C_f$')

# FDS Data
h_leg = []
for q in range(lnx):
   df = M[q]

   # Column names assumed from MATLAB script: 'u_tau-x' and 'u_tau'
   x = df['u_tau-x'].values
   I = x > 0
   x = x[I]

   u_tau = df['u_tau'].values[I]
   # tau_w = rho * u_tau**2 # Not used for Cf calculation here

   # Column name assumed from MATLAB script: 'u_wall'
   u_wall = df['u_wall'].values[I]

   dx_q = h / nx[q]
   Cf_fds = visc * u_wall / (0.5 * dx_q) / (0.5 * rho * U_0**2)

   fdsplotlib.plot_to_fig(x_data=x/h, y_data= Cf_fds, marker_style=fds_marker[q],data_markevery=int(len(Cf_fds)/(14+q)),
         figure_handle=fig,
         data_label=fds_key[q])

plotname = pltdir + 'backward_facing_step_Cf.pdf'
plt.savefig(plotname, format='pdf')
plt.close()

# --- 2. Pressure coefficient along bottom of channel (Cp) ---

# Exp Data (Cp-x/h)
xoh_cp = D['Cp-x/h'].values
Cp_exp = D['Cp'].values

# Filter data where x/h <= 20
I = xoh_cp <= 20
xoh_cp = xoh_cp[I]
Cp_exp = Cp_exp[I]

fig = fdsplotlib.plot_to_fig(x_data=xoh_cp, y_data=Cp_exp, marker_style='ko',
      revision_label=version_string,x_min=0,x_max=20,y_min=-0.1,y_max=0.25,
      data_label='J&D',
      x_label=r'$x/h$',
      y_label=r'$C_p$')

# FDS Data
h_leg2 = []
for q in range(lnx):
   df = M[q]

   x = df['cp-x'].values
   I = x > 0
   x = x[I]

   Cp_fds_raw = df['cp'].values[I]

   # Cp_fds = Cp / U_0^2
   Cp_fds = Cp_fds_raw / (U_0**2)

   # Normalize Cp data to match experimental end point
   # Cp_fds = Cp_fds + Cp(end) - Cp_fds(end);
   Cp_fds = Cp_fds + Cp_exp[-1] - Cp_fds[-1]


   fdsplotlib.plot_to_fig(x_data=x/h, y_data= Cp_fds, marker_style=fds_marker[q],data_markevery=int(len(Cp_fds)/(14+q)),
         figure_handle=fig,
         data_label=fds_key[q])

plotname = pltdir + 'backward_facing_step_Cp.pdf'
plt.savefig(plotname, format='pdf')
plt.close()

# --- 3. Wall-normal streamwise velocity profiles (U) ---

fig, sp1 = tight_subplot(1, 4, gap=(0.01, 0.01), marg_h=(0.142, 0.055), marg_w=(0.108, 0.01))
fig.suptitle(version_string,x=0.8,y=0.98,fontsize=14)

for i in range(4): # Loop over x_loc: -3, 4, 6, 10
   ax = sp1[i]
   x_l = x_loc[i]

   # Exp Data
   z_col = f'z {x_l}'
   u_col = f'U {x_l}'

   z_exp = D[z_col].values
   u_data = D[u_col].values

   I = z_exp > 0
   z_exp = z_exp[I] * 0.001 # Convert mm to m
   u_data = u_data[I]

   if i == 0: # x/h = -3
      z_exp += h # Add step height to z coordinate

   ax.plot(u_data / U_0, z_exp / h, 'ko', markersize=10)

   # FDS Data
   for q in range(lnx):
      df = M[q]

      z_fds_col = f'U-VEL {x_l}-z'
      u_fds_col = f'U-VEL {x_l}'

      z_fds = df[z_fds_col].values
      u_vel = df[u_fds_col].values

      I = z_fds > 0
      z_fds = z_fds[I]
      u_vel = u_vel[I]

      # Plot line and markers
      ax.plot(
         u_vel/ U_0,
         z_fds/ h,
         marker=subsym[q],
         markevery=int(len(z_fds) /(14+q)),
         linestyle='-',
         color=subcol[q],
         linewidth=1,
         markersize=10
      )

   ax.set_xlim(-0.2, 1.0)
   ax.set_ylim(0, 3.5)
   ax.tick_params(axis='both', which='major', labelsize=14)
   ax.text(0.01, 3.2, f'$x/h$={x_l}', transform=ax.transData, fontsize=14)

   if i == 0:
      ax.set_ylabel(r'$z/h$', fontsize=16)
      ax.set_xlabel(r'$\langle u \rangle/U_0$', fontsize=16)
      # Add legend only to the first plot in this set
      legend_labels = ['J&D']
      for k in fds_key:
         legend_labels.append(k)
      ax.legend(legend_labels, loc='best', frameon=False, fontsize=10)
   else:
      ax.set_yticklabels([])
      ax.set_xticklabels([])

fig.subplots_adjust(left=0.108, right=1.0 - 0.01, bottom=0.055, top=1.0 - 0.142)
fig.savefig(os.path.join(pltdir, 'backward_facing_step_U.pdf'))
plt.close(fig)

# --- 4. Wall-normal vertical velocity profiles (W) ---

fig, sp1 = tight_subplot(1, 4, gap=(0.01, 0.01), marg_h=(0.142, 0.055), marg_w=(0.108, 0.01))
fig.suptitle(version_string,x=0.8,y=0.98,fontsize=14)

for i in range(4):
   ax = sp1[i]
   x_l = x_loc[i]

   # Exp Data
   z_col = f'z {x_l}'
   w_col = f'W {x_l}'

   z_exp = D[z_col].values
   w_data = D[w_col].values

   I = z_exp > 0
   z_exp = z_exp[I] * 0.001
   w_data = w_data[I]

   if i == 0:
      z_exp += h

   ax.plot(w_data / U_0, z_exp / h, 'ko', markersize=10)

   # FDS Data
   for q in range(lnx):
      df = M[q]

      z_fds_col = f'W-VEL {x_l}-z'
      w_fds_col = f'W-VEL {x_l}'

      z_fds = df[z_fds_col].values
      w_vel = df[w_fds_col].values

      I = z_fds > 0
      z_fds = z_fds[I]
      w_vel = w_vel[I]

      ax.plot(
         w_vel/ U_0,
         z_fds/ h,
         marker=subsym[q],
         markevery=int(len(z_fds) /(14+q)),
         linestyle='-',
         color=subcol[q],
         linewidth=1,
         markersize=10
      )

   ax.set_xlim(-0.1, 0.1)
   ax.set_ylim(0, 3.5)
   ax.tick_params(axis='both', which='major', labelsize=14)
   ax.text(0.01, 3.2, f'$x/h$={x_l}', transform=ax.transData, fontsize=14)

   if i == 0:
      ax.set_ylabel(r'$z/h$', fontsize=16)
      ax.set_xlabel(r'$\langle w \rangle/U_0$', fontsize=16)
      # Add legend only to the first plot in this set
      legend_labels = ['J&D']
      for k in fds_key:
         legend_labels.append(k)
      ax.legend(legend_labels, loc='best', frameon=False, fontsize=10)
   else:
      ax.set_yticklabels([])
      ax.set_xticklabels([])

fig.subplots_adjust(left=0.108, right=1.0 - 0.01, bottom=0.055, top=1.0 - 0.142)
fig.savefig(os.path.join(pltdir, 'backward_facing_step_W.pdf'))
plt.close(fig)

# --- 5. Wall-normal uu profiles ($\langle uu \rangle$) ---

fig, sp1 = tight_subplot(1, 4, gap=(0.01, 0.01), marg_h=(0.142, 0.055), marg_w=(0.108, 0.01))
fig.suptitle(version_string,x=0.8,y=0.98,fontsize=14)

for i in range(4):
   ax = sp1[i]
   x_l = x_loc[i]

   # Exp Data
   z_col = f'z {x_l}'
   uu_col = f'uu {x_l}'

   z_exp = D[z_col].values
   uu_data = D[uu_col].values

   I = z_exp > 0
   z_exp = z_exp[I] * 0.001
   uu_data = uu_data[I]

   if i == 0:
      z_exp += h

   # No error bars plotted in the original (herrorbar was commented out)
   ax.plot(uu_data / U_0**2, z_exp / h, 'ko', markersize=10)

   # FDS Data
   for q in range(lnx):
      df = M[q]

      z_fds_col = f'uu {x_l}-z'
      uu_fds_col = f'uu {x_l}'

      z_fds = df[z_fds_col].values
      uu_fds = df[uu_fds_col].values

      I = z_fds > 0
      z_fds = z_fds[I]
      uu_fds = uu_fds[I]

      ax.plot(
         uu_fds/ U_0**2,
         z_fds/ h,
         marker=subsym[q],
         markevery=int(len(z_fds) /(14+q)),
         linestyle='-',
         color=subcol[q],
         linewidth=1,
         markersize=10
      )

   ax.set_xlim(0, 0.04)
   ax.set_ylim(0, 3.5)
   ax.tick_params(axis='both', which='major', labelsize=14)
   ax.text(0.01, 3.2, f'$x/h={x_l}$', transform=ax.transData, fontsize=14)

   if i == 0:
      ax.set_ylabel(r'$z/h$', fontsize=16)
      ax.set_xlabel(r'$\langle uu \rangle/U_0^2$', fontsize=16)
      # Add legend only to the first plot in this set
      legend_labels = ['J&D']
      for k in fds_key:
         legend_labels.append(k)
      ax.legend(legend_labels, loc='best', frameon=False, fontsize=10)
   else:
      ax.set_yticklabels([])
      ax.set_xticklabels([])

fig.subplots_adjust(left=0.108, right=1.0 - 0.01, bottom=0.055, top=1.0 - 0.142)
fig.savefig(os.path.join(pltdir, 'backward_facing_step_uu.pdf'))
plt.close(fig)

# --- 6. Wall-normal ww profiles ($\langle ww \rangle$) ---

fig, sp1 = tight_subplot(1, 4, gap=(0.01, 0.01), marg_h=(0.142, 0.055), marg_w=(0.108, 0.01))
fig.suptitle(version_string,x=0.8,y=0.98,fontsize=14)

for i in range(4):
   ax = sp1[i]
   x_l = x_loc[i]

   # Exp Data
   z_col = f'z {x_l}'
   ww_col = f'ww {x_l}'

   z_exp = D[z_col].values
   ww_data = D[ww_col].values

   I = z_exp > 0
   z_exp = z_exp[I] * 0.001
   ww_data = ww_data[I]

   if i == 0:
      z_exp += h

   ax.plot(ww_data / U_0**2, z_exp / h, 'ko', markersize=10)

   # FDS Data
   for q in range(lnx):
      df = M[q]

      z_fds_col = f'ww {x_l}-z'
      ww_fds_col = f'ww {x_l}'

      z_fds = df[z_fds_col].values
      ww_fds = df[ww_fds_col].values

      I = z_fds > 0
      z_fds = z_fds[I]
      ww_fds = ww_fds[I]

      ax.plot(
         ww_fds/ U_0**2,
         z_fds/ h,
         marker=subsym[q],
         markevery=int(len(z_fds) /(14+q)),
         linestyle='-',
         color=subcol[q],
         linewidth=1,
         markersize=10
      )

   ax.set_xlim(0, 0.04)
   ax.set_ylim(0, 3.5)
   ax.tick_params(axis='both', which='major', labelsize=14)
   ax.text(0.01, 3.2, f'$x/h$={x_l}', transform=ax.transData, fontsize=14)

   if i == 0:
      ax.set_ylabel(r'$z/h$', fontsize=16)
      ax.set_xlabel(r'$\langle ww \rangle/U_0^2$', fontsize=16)
      # Add legend only to the first plot in this set
      legend_labels = ['J&D']
      for k in fds_key:
         legend_labels.append(k)
      ax.legend(legend_labels, loc='best', frameon=False, fontsize=10)
   else:
      ax.set_yticklabels([])
      ax.set_xticklabels([])

fig.subplots_adjust(left=0.108, right=1.0 - 0.01, bottom=0.055, top=1.0 - 0.142)
fig.savefig(os.path.join(pltdir, 'backward_facing_step_ww.pdf'))
plt.close(fig)

fig, sp1 = tight_subplot(1, 4, gap=(0.01, 0.01), marg_h=(0.142, 0.055), marg_w=(0.108, 0.01))
fig.suptitle(version_string,x=0.8,y=0.98,fontsize=14)

# --- 7. Wall-normal uw profiles ($\langle uw \rangle$) ---

for i in range(4):
   ax = sp1[i]
   x_l = x_loc[i]

   # Exp Data
   z_col = f'z {x_l}'
   uw_col = f'uw {x_l}'

   z_exp = D[z_col].values
   uw_data = D[uw_col].values

   I = z_exp > 0
   z_exp = z_exp[I] * 0.001
   uw_data = uw_data[I]

   if i == 0:
      z_exp += h

   ax.plot(uw_data / U_0**2, z_exp / h, 'ko', markersize=10)

   # FDS Data
   h_leg = []
   for q in range(lnx):
      df = M[q]

      z_fds_col = f'uw {x_l}-z'
      uw_fds_col = f'uw {x_l}'

      z_fds = df[z_fds_col].values
      uw_fds_raw = df[uw_fds_col].values

      I = z_fds > 0
      z_fds = z_fds[I]
      uw_fds_raw = uw_fds_raw[I]

      # Original MATLAB applies -1 multiplication: uw_fds = -1*M{q}.data(I,j);
      uw_fds = -1 * uw_fds_raw
      ax.plot(
         uw_fds/ U_0**2,
         z_fds/ h,
         marker=subsym[q],
         markevery=int(len(z_fds) /(14+q)),
         linestyle='-',
         color=subcol[q],
         linewidth=1,
         markersize=10
      )

   ax.set_xlim(0, 0.04) # Note: Original MATLAB had axis([0 .04 0 3.5])
   ax.set_ylim(0, 3.5)
   ax.tick_params(axis='both', which='major', labelsize=14)
   ax.text(0.01, 3.2, f'$x/h$={x_l}', transform=ax.transData, fontsize=14)

   if i == 0:
      ax.set_ylabel(r'$z/h$', fontsize=16)
      ax.set_xlabel(r'$\langle uw \rangle/U_0^2$', fontsize=16)
      # Add legend only to the first plot in this set
      legend_labels = ['J&D']
      for k in fds_key:
         legend_labels.append(k)
      ax.legend(legend_labels, loc='best', frameon=False, fontsize=10)
   else:
      ax.set_yticklabels([])
      ax.set_xticklabels([])

fig.subplots_adjust(left=0.108, right=1.0 - 0.01, bottom=0.055, top=1.0 - 0.142)
fig.savefig(os.path.join(pltdir, 'backward_facing_step_uw.pdf'))
plt.close(fig)

