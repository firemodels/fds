# McDermott and Henderson
# 26 Oct 2017 (updated 11-27-2022 RJM)
#% natconh.m
#
# Converted by Floyd
# 10/21/2025


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import fdsplotlib
import os
import sys

results_dir = '../../../out/Convection/'
casename = [
'natconh_1',
'natconh_2',
'natconh_3',
'natconh_4',
'natconh_5',
'natconh_6',
'natconh_7',
'natconh_8',
'natconh_9',
'natconh_10',
'natconh_11',
'natconh_12',
'natconh_13',
'natconh_14',
'natconh_15',
'natconh_16',
'natconh_17',
'natconh_18',
'natconh_19'
]

res = ['8', '16', '32']
res_cg = ['8', '16']

# calculations below were used for input file setup

g = 9.80665
S = np.array([.002, .002, .02, .02, .02, .02, .02, .2, .2, .2, .2, 2, 2, 2, 2, 20, 20, 20, 20])
T1 = np.array([303.15, 600.15, 295.15, 297.15, 303.15, 333.15, 503.15, 295.15, 303.15, 333.15, 503.15, 295.15, 303.15, 333.15, 503.15, 295.15, 303.15, 333.15, 503.15])
T2 = 293.15
Tm = 0.5 * (T1 + T2)
beta = 1.0 / Tm
MW = 28.85476 # FDS 'LJ AIR'
P0 = 101325
rho = P0 * MW / (8341.472 * Tm)
mu = 1.8216e-5
cp = 1000
k = 0.018216 # for Pr=1 fluid

Pr = cp * mu / k
nu = mu / rho
alpha = k / (rho * cp)

Ra = (g * beta * (T1 - T2) * S**3) / (alpha * nu)

# see J.P. Holman p. 361 for correlations

Ra_Limit_1 = 1183
Ra_Limit_2 = 5049
Ra_Limit_3 = 3.10508e6

# Check for files

skip_case = False
for i in range(len(casename)):
   for j in range(len(res)):
      devc_filepath = os.path.join(results_dir, f"{casename[i]}_{res[j]}_devc.csv")
      if not os.path.exists(devc_filepath):
         skip_case = True
         print('Error: File ', devc_filepath, ' does not exist. Skipping case.')
      line_filepath = os.path.join(results_dir, f"{casename[i]}_{res[j]}_line.csv")
      if not os.path.exists(line_filepath):
         skip_case = True
         print('Error: File ', line_filepath, ' does not exist. Skipping case.')
   for j in range(len(res_cg)):
      devc_filepath = os.path.join(results_dir, f"{casename[i]}_{res_cg[j]}_rot_18_devc.csv")
      if not os.path.exists(devc_filepath):
         skip_case = True
         print('Error: File ', devc_filepath, ' does not exist. Skipping case.')

if skip_case: quit()

r0 = []
r1 = []
r2 = []
r3 = []
Nu = np.zeros_like(Ra)

for i in range(len(Ra)):
   if Ra[i] < Ra_Limit_1:
      r0.append(i)
      Nu[i] = 1
   elif Ra[i] > Ra_Limit_1 and Ra[i] <= Ra_Limit_2:
      r1.append(i)
      Nu[i] = 0.059 * Ra[i]**(0.4)
   elif Ra[i] > Ra_Limit_2 and Ra[i] <= Ra_Limit_3:
      r2.append(i)
      Nu[i] = 0.212 * Ra[i]**(0.25)
   elif Ra[i] > Ra_Limit_3:
      r3.append(i)
      Nu[i] = 0.061 * Ra[i]**(1/3)

RAYLEIGH = np.logspace(0, 14, 1000)
NUSSELT = np.zeros_like(RAYLEIGH)
for i in range(len(RAYLEIGH)):
   if RAYLEIGH[i] < Ra_Limit_1:
      NUSSELT[i] = 1
   elif RAYLEIGH[i] > Ra_Limit_1 and RAYLEIGH[i] <= Ra_Limit_2:
      NUSSELT[i] = 0.059 * RAYLEIGH[i]**(0.4)
   elif RAYLEIGH[i] > Ra_Limit_2 and RAYLEIGH[i] <= Ra_Limit_3:
      NUSSELT[i] = 0.212 * RAYLEIGH[i]**(0.25)
   elif RAYLEIGH[i] > Ra_Limit_3:
      NUSSELT[i] = 0.061 * RAYLEIGH[i]**(1/3)

# --- Plot 1: Correlation and FDS Results ---

Git_Filename = os.path.join(results_dir, 'natconh_19_32_git.txt')
version_string = fdsplotlib.get_version_string(Git_Filename)

fig = fdsplotlib.plot_to_fig(x_data=RAYLEIGH, y_data=NUSSELT, marker_style='k-',plot_type='loglog',
      revision_label=version_string,x_min=1,x_max=1e15,y_min=0.5,y_max=1e4,
      data_label=r'Correlation, Nu = $C$ Ra$^n$',
      x_label=r'Rayleigh Number',
      y_label=r'Nusselt Number')

delta = S
T = T1

marker_style = ['r^', 'gs', 'b+']
qcolhdrs = ['Q1-1', 'Q1-2', 'Q1-3', 'Q1-4', 'Q1-5', 'Q2-1', 'Q2-2', 'Q2-3', 'Q2-4', 'Q2-5']

for j in range(len(res)):
   Ra_plot = np.zeros(len(Ra))
   Nu_plot = np.zeros(len(Ra))
   for i in range(len(delta)):

      devc_filepath = os.path.join(results_dir, f"{casename[i]}_{res[j]}_devc.csv")
      M_devc = pd.read_csv(devc_filepath, skiprows=[1])

      t = M_devc.iloc[:, 0].values
      Q1 = M_devc.iloc[:, 1].values
      Q2 = M_devc.iloc[:, 2].values # Third column
      rho = M_devc.iloc[-1, 4]
      alpha_calc = k / (rho * cp)
      nu_calc = mu / rho
      b = 2.0 / (T[i] + T2)
      Ra_plot[i] = (g * b * (T[i] - T2) * delta[i]**3) / (alpha_calc * nu_calc)

      line_filepath = os.path.join(results_dir, f"{casename[i]}_{res[j]}_line.csv")
      M_line = pd.read_csv(line_filepath, skiprows=1)

      cols = [M_line.columns.get_loc(hdr) for hdr in qcolhdrs if hdr in M_line.columns]

      Q = np.mean(np.mean(np.abs(M_line.iloc[:, cols].values) * 1000)) # heat flow, W/m2

      Nu_plot[i] = Q * (delta[i] / k) / (T[i] - T2)

   plotlabel = fr'FDS $H/\Delta x$={res[j]}'
   fdsplotlib.plot_to_fig(x_data=Ra_plot, y_data=Nu_plot, marker_style=marker_style[j],marker_fill_color='none',
      figure_handle=fig,
         data_label = plotlabel)

plt.savefig('../../Manuals/FDS_Verification_Guide/SCRIPT_FIGURES/natural_convection_horizontal_enclosure.pdf', format='pdf')

# --- Plot 2: Complex geometry cases (rotated grid 18 degrees) ---
Git_Filename = os.path.join(results_dir, 'natconh_19_16_rot_18_git.txt')
version_string = fdsplotlib.get_version_string(Git_Filename)

fig = fdsplotlib.plot_to_fig(x_data=RAYLEIGH, y_data=NUSSELT, marker_style='k-',plot_type='loglog',
      revision_label=version_string,x_min=1,x_max=1e15,y_min=0.5,y_max=1e4,
      data_label=r'Correlation, Nu = C Ra$^n$',
      x_label=r'Rayleigh Number',
      y_label=r'Nusselt Number')

for j in range(len(res_cg)):
   for i in range(len(delta)):

      devc_filepath = os.path.join(results_dir, f"{casename[i]}_{res[j]}_rot_18_devc.csv")
      M_devc_geom = pd.read_csv(devc_filepath, skiprows=1)

      t = M_devc_geom.iloc[:, 0].values # Time column
      Q1 = M_devc_geom['Q1-0'].values
      Q2 = M_devc_geom['Q2-0'].values
      rho = M_devc_geom['rho'].iloc[-1]

      alpha_calc = k / (rho * cp)
      nu_calc = mu / rho
      b = 2.0 / (T[i] + T2)
      Ra_plot[i] = (g * b * (T[i] - T2) * delta[i]**3) / (alpha_calc * nu_calc)

      qrange = np.where(t > t[-1] / 2)[0]
      A = (delta[i] * 8)**2
      q = np.mean(np.abs(Q2[qrange])) * 1000 / A # heat flux, W/m2

      Nu_plot[i] = q * (delta[i] / k) / (T[i] - T2)

   plotlabel = fr'GEOM $H/\Delta x$={res_cg[j]}'
   fdsplotlib.plot_to_fig(x_data=Ra_plot, y_data=Nu_plot, marker_style=marker_style[j],marker_fill_color='none',
      figure_handle=fig,
         data_label = plotlabel)

plt.savefig('../../Manuals/FDS_Verification_Guide/SCRIPT_FIGURES/natconh_geom.pdf', format='pdf')
