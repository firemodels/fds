# McDermott
# 10-7-14
# hot_layer_collapse.m
#
# H. Baum. "Collapse of a Hot Layer in a Micro-Gravity Environment", Sep 2, 2014 (personal notes)
#
# Converted by Floyd
# 10/16/2025

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import fdsplotlib
from scipy.special import erfc, erfinv

# --- Configuration (Approximation of MATLAB external variables) ---
datadir = '../../Verification/Flowfields/'
plotdir = '../../Manuals/FDS_Verification_Guide/SCRIPT_FIGURES/'
filename = 'hot_layer_360_devc.csv'
git_file = datadir+ 'hot_layer_360_git.txt'

skip_case = False

name = datadir+filename
if not os.path.exists(name):
   skip_case = True
   print('Error: File ', filename, ' does not exist. Skipping case.')

if skip_case: quit()

error_tolerance = 0.01

L = 2.0      # height of domain 2*d where d is the layer depth at t=0
N = 360      # number of cells in vertical direction
T_0 = 293.15 # cold wall temperature (K)
T_h = 1172.6 # hot layer initial temperature (K)
k_0 = 1.0    # W/m/K
rho_0 = 1.0  # kg/m^3
Cp = 1000.0  # J/kg/K
d = 1.0      # m
v0 = k_0 / (rho_0 * Cp * d) # m/s

dz = L / N

t = np.array([.2, .4, .6, .8, 1, 2, 4, 6, 8, 10])
tau = t * 1e-3
marker_style = ['ko', 'mo', 'ro', 'go', 'bo']
exact_soln_style = ['k-', 'm-', 'r-', 'g-', 'b-']
legend_entries = [r'$\tau = 2 \times 10^{-4}$', r'$\tau = 4 \times 10^{-4}$', r'$\tau = 6 \times 10^{-4}$', r'$\tau = 8 \times 10^{-4}$', r'$\tau = 10 \times 10^{-4}$']
legend_entries_2 = [r'$\tau = 2 \times 10^{-3}$', r'$\tau = 4 \times 10^{-3}$', r'$\tau = 6 \times 10^{-3}$', r'$\tau = 8 \times 10^{-3}$', r'$\tau = 10 \times 10^{-3}$']

def f(lambda_in, tau_in, a_in):
   term1 = -2 * np.sqrt(tau_in / np.pi) * np.exp(-(a_in - lambda_in)**2 / (4 * tau_in))
   term2 = (a_in + lambda_in) * erfc((a_in + lambda_in) / (2 * np.sqrt(tau_in)))
   return term1 + term2

def F(lambda_in, tau_in, a_in):
   return f(lambda_in, tau_in, a_in) - f(0, tau_in, a_in)

def U(lambda_in, tau_in, T_h_in, T_0_in):
   arg1 = lambda_in / (2 * np.sqrt(tau_in))
   arg2 = (lambda_in - T_0_in / T_h_in) / (2 * np.sqrt(tau_in))
   arg3 = (lambda_in + T_0_in / T_h_in) / (2 * np.sqrt(tau_in))

   return -erfc(arg1) + 0.5 * (erfc(arg2) + erfc(arg3))

def Theta(lambda_in, tau_in, T_h_in, T_0_in):
   return 1.0 + ((T_h_in - T_0_in) / T_0_in) * U(lambda_in, tau_in, T_h_in, T_0_in)

def G(lambda_in, tau_in, a_in):
   term1 = -np.exp(-(a_in + lambda_in)**2 / (4 * tau_in))
   term2 = np.exp(-a_in**2 / (4 * tau_in))
   return (1 / np.sqrt(np.pi * tau_in)) * (term1 + term2)

def V(lambda_in, tau_in, T_h_in, T_0_in):
   term1 = -G(lambda_in, tau_in, 0)
   term2 = 0.5 * (G(lambda_in, tau_in, T_0_in / T_h_in) + G(lambda_in, tau_in, -T_0_in / T_h_in))
   return ((T_h_in - T_0_in) / T_0_in) * (term1 + term2)

def Y(lambda_in, tau_in, T_h_in, T_0_in):
   term1 = -F(lambda_in, tau_in, 0)
   term2 = 0.5 * (F(lambda_in, tau_in, T_0_in / T_h_in) + F(lambda_in, tau_in, -T_0_in / T_h_in))
   return lambda_in + ((T_h_in - T_0_in) / T_0_in) * (term1 + term2)


lambda_range = np.linspace(0, 10 * T_0 / T_h, 1000)

df = pd.read_csv(datadir+filename, skiprows=2, header=None)

# Find row indices 'J' corresponding to times 't'
J = []
for time_t in t:
   try:
      index_j = df[df[0] >= time_t].index[0]
      J.append(index_j)
   except IndexError:
      print(f"Warning: Time {time_t} not found in data. Stopping.")
      break
if not J:
   exit()

T_range_indices = slice(1, 41)
W_range_indices = slice(41, 81)

dz_40 = L / 40.0
dz_360 = L / 360.0
z_T = np.arange(dz_40 / 2, L, dz_40)
z_W = np.arange(dz_40 / 2 + dz_360 / 2, L, dz_40)

if len(z_T) > 40:
   z_T = z_T[:40]
if len(z_W) > 40:
   z_W = z_W[:40]

version_string = fdsplotlib.get_version_string(git_file)

for j in range(5): # j=0 to 4 (corresponding to tau[0] to tau[4])
   T_C = df.iloc[J[j], T_range_indices].values
   T_nondim = (T_C + 273.15) / T_0

   if j==0:
      fig = fdsplotlib.plot_to_fig(x_data=T_nondim, y_data=z_T, marker_style=marker_style[j],
      revision_label=version_string,x_min=1,x_max=4,y_min=0,y_max=1.2,
      data_label=legend_entries[j],
      legend_location='center left',
      y_label=r'$T/T_0$',
      x_label=r'$y/d$')
   else:
      fdsplotlib.plot_to_fig(x_data=T_nondim, y_data=z_T, marker_style=marker_style[j],
         figure_handle=fig,
         data_label=legend_entries[j],)

   fdsplotlib.plot_to_fig(x_data=Theta(lambda_range, tau[j], T_h, T_0), y_data= Y(lambda_range, tau[j], T_h, T_0), marker_style=exact_soln_style[j],
      figure_handle=fig)

plotfile = plotdir + 'hot_layer_temp_1.pdf'
plt.savefig(plotfile, format='pdf')
plt.close()

# --- Compute Temperature Error at tau[4] (j=5 in MATLAB, index 4 in Python) ---
j_err = 4
T_C_err = df.iloc[J[j_err], T_range_indices].values
T_nondim_err = (T_C_err + 273.15) / T_0

# map nondimensional position: Find index 'I' where Y >= z_T
I = []
Y_analytic = Y(lambda_range, tau[j_err], T_h, T_0)
for z_val in z_T:
   # Find the first index in lambda_ where Y_analytic >= z_val
   try:
      I.append(np.where(Y_analytic >= z_val)[0][0])
   except IndexError:
      print(f"Warning: Analytical Y does not cover z_T={z_val}")
      I.append(I[-1] if I else 0)

# Ensure I has the same length as z_T
I = np.array(I[:len(z_T)])

T_analytic_at_FDS_pos = Theta(lambda_range[I], tau[j_err], T_h, T_0)

Error_T1 = np.linalg.norm(T_nondim_err - T_analytic_at_FDS_pos) / (np.max(T_nondim_err) * len(T_nondim_err))

if Error_T1 > error_tolerance:
    print(f"Python Warning: hot_layer_360.fds Temp_1 Error = {Error_T1:.4e}")

# --- Plot 2: Velocity v/v_0 (Short Times) ---
fig2, ax2 = plt.subplots()

for j in range(5):
   # Extract FDS Velocity, convert to v/v_0: W / v0
   W = df.iloc[J[j], W_range_indices].values
   W_nondim = W / v0

   if j==0:
      fig = fdsplotlib.plot_to_fig(x_data=W_nondim, y_data=z_W, marker_style=marker_style[j],
      revision_label=version_string,x_min=-200,x_max=0,y_min=0,y_max=1.2,
      data_label=legend_entries[j],
      y_label=r'$v/v_0$',
      x_label=r'$y/d$')
   else:
      fdsplotlib.plot_to_fig(x_data=W_nondim, y_data=z_W, marker_style=marker_style[j],
         figure_handle=fig,
         data_label=legend_entries[j],)

   fdsplotlib.plot_to_fig(x_data=V(lambda_range, tau[j], T_h, T_0), y_data= Y(lambda_range, tau[j], T_h, T_0), marker_style=exact_soln_style[j],
      figure_handle=fig)

plotfile = plotdir + 'hot_layer_vel_1.pdf'
plt.savefig(plotfile, format='pdf')
plt.close()

# --- Compute Velocity Error at tau[4] (j=5) ---
j_err = 4
W_err = df.iloc[J[j_err], W_range_indices].values
W_nondim_err = W_err / v0

# map nondimensional position: Find index 'I' where Y >= z_W
I = []
Y_analytic = Y(lambda_range, tau[j_err], T_h, T_0)
for z_val in z_W:
   try:
      I.append(np.where(Y_analytic >= z_val)[0][0])
   except IndexError:
      I.append(I[-1] if I else 0)

I = np.array(I[:len(z_W)])

V_analytic_at_FDS_pos = V(lambda_range[I], tau[j_err], T_h, T_0)

Error_V1 = np.linalg.norm(W_nondim_err - V_analytic_at_FDS_pos) / (np.max(np.abs(W_nondim_err)) * len(W_nondim_err))

if Error_V1 > error_tolerance:
    print(f"Python Warning: hot_layer_360.fds Vel_1 Error = {Error_V1:.4e}")

# ----------------------------------------------------------------------
# --- Second set of Plots (Longer Times: tau[5] to tau[9]) ---
# ----------------------------------------------------------------------

# --- Plot 3: Temperature T/T_0 (Longer Times) ---
fig3, ax3 = plt.subplots()

for j in range(5, 10): # j=5 to 9 (corresponding to tau[5] to tau[9])
   jj = j - 5 # new index for markers/legend (0 to 4)

   T_C = df.iloc[J[j], T_range_indices].values
   T_nondim = (T_C + 273.15) / T_0

   if j==5:
      fig = fdsplotlib.plot_to_fig(x_data=T_nondim, y_data=z_T, marker_style=marker_style[jj],
      revision_label=version_string,x_min=1,x_max=4,y_min=0,y_max=1.2,
      data_label=legend_entries[jj],
      y_label=r'$T/T_0$',
      x_label=r'$y/d$')
   else:
      fdsplotlib.plot_to_fig(x_data=T_nondim, y_data=z_T, marker_style=marker_style[jj],
         figure_handle=fig,
         data_label=legend_entries[jj],)

   fdsplotlib.plot_to_fig(x_data=Theta(lambda_range, tau[j], T_h, T_0), y_data= Y(lambda_range, tau[j], T_h, T_0), marker_style=exact_soln_style[jj],
      figure_handle=fig)

plotfile = plotdir + 'hot_layer_temp_2.pdf'
plt.savefig(plotfile, format='pdf')
plt.close()

# --- Compute Temperature Error at tau[9] (j=10) ---
j_err = 9
T_C_err = df.iloc[J[j_err], T_range_indices].values
T_nondim_err = (T_C_err + 273.15) / T_0

# map nondimensional position
I = []
Y_analytic = Y(lambda_range, tau[j_err], T_h, T_0)
for z_val in z_T:
   try:
      I.append(np.where(Y_analytic >= z_val)[0][0])
   except IndexError:
      I.append(I[-1] if I else 0)

I = np.array(I[:len(z_T)])

T_analytic_at_FDS_pos = Theta(lambda_range[I], tau[j_err], T_h, T_0)
Error_T2 = np.linalg.norm(T_nondim_err - T_analytic_at_FDS_pos) / (np.max(T_nondim_err) * len(T_nondim_err))

if Error_T2 > error_tolerance:
    print(f"Python Warning: hot_layer_360.fds Temp_2 Error = {Error_T2:.4e}")

# --- Plot 4: Velocity v/v_0 (Longer Times) ---
fig4, ax4 = plt.subplots()

for j in range(5, 10):
   jj = j - 5

   W = df.iloc[J[j], W_range_indices].values
   W_nondim = W / v0

   if j==5:
      fig = fdsplotlib.plot_to_fig(x_data=W_nondim, y_data=z_W, marker_style=marker_style[jj],
      revision_label=version_string,x_min=-60,x_max=0,y_min=0,y_max=1.6,
      data_label=legend_entries[jj],
      y_label=r'$v/v_0$',
      x_label=r'$y/d$')
   else:
      fdsplotlib.plot_to_fig(x_data=W_nondim, y_data=z_W, marker_style=marker_style[jj],
         figure_handle=fig,
         data_label=legend_entries[jj],)

   fdsplotlib.plot_to_fig(x_data=V(lambda_range, tau[j], T_h, T_0), y_data= Y(lambda_range, tau[j], T_h, T_0), marker_style=exact_soln_style[jj],
      figure_handle=fig)

plotfile = plotdir + 'hot_layer_vel_2.pdf'
plt.savefig(plotfile, format='pdf')
plt.close()

# --- Compute Velocity Error at tau[9] (j=10) ---
j_err = 9
W_err = df.iloc[J[j_err], W_range_indices].values
W_nondim_err = W_err / v0

# map nondimensional position
I = []
Y_analytic = Y(lambda_range, tau[j_err], T_h, T_0)
for z_val in z_W:
   try:
      I.append(np.where(Y_analytic >= z_val)[0][0])
   except IndexError:
      I.append(I[-1] if I else 0)

I = np.array(I[:len(z_W)])

V_analytic_at_FDS_pos = V(lambda_range[I], tau[j_err], T_h, T_0)
Error_V2 = np.linalg.norm(W_nondim_err - V_analytic_at_FDS_pos) / (np.max(np.abs(W_nondim_err)) * len(W_nondim_err))

if Error_V2 > error_tolerance:
   print(f"Python Warning: hot_layer_360.fds Vel_2 Error = {Error_V2:.4e}")