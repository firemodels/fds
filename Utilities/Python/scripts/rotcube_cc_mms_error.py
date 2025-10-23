# Vanella rotcube_cc_mms_error.m
# 05-24-2019
# Computes L1, L2, and Linf errors for the Rotated Cube MMS case.
#
# Converted by Floyd
# 10/17/2025

import numpy as np
import pandas as pd
import os
import math
import matplotlib.pyplot as plt
import fdsplotlib


mu = 0.01
D  = 0.01
rho = 1.00
L = np.pi # Cube side length.

nwave = 1. # Wave number on analytical solution.
gam = np.pi/2
Az = 0.1
meanz = 0.15
displ = np.pi
dispxy = np.array([-1/2 * displ, -1/2 * displ])

datadir = '../../Verification/Complex_Geometry/'
plotdir = '../../Manuals/FDS_Verification_Guide/SCRIPT_FIGURES/'
gitname = 'rotated_cube_45deg_256_stm_git.txt'
git_file = datadir+gitname
version_string = fdsplotlib.get_version_string(git_file)

def ann_z(x, y, t):
   return (
      Az / 3. * np.sin(t) * (1. - np.cos(2. * nwave * (x - gam))) * (1. - np.cos(2. * nwave * (y - gam)))
      - Az / 3. * np.sin(t)
      + meanz
   )

def ann_u(x, y, t):
   return -np.sin(t) * np.sin(nwave * x)**2 * np.sin(2. * nwave * y)

def ann_v(x, y, t):
   return np.sin(t) * np.sin(2. * nwave * x) * np.sin(nwave * y)**2

def ann_p(x, y, t):
   return np.sin(t) / 4. * (2. + np.cos(2. * nwave * x)) * (2. + np.cos(2. * nwave * y)) - np.sin(t)

def ann_H(x, y, t):
   u = ann_u(x, y, t)
   v = ann_v(x, y, t)
   return (u**2 + v**2) / 2 + ann_p(x, y, t) / rho


# Function to read the MATLAB-style CSV file header
def read_matlab_csv(filepath):
   try:
      # Read the first line (metadata)
      with open(filepath, 'r') as f:
         header_str = f.readline().strip()

      # MATLAB's str2num is flexible; we try to handle space/comma separation
      # and ensure we get 6 numbers.
      vec = [float(x) for x in header_str.replace(',', ' ').split() if x.strip()]

      if len(vec) < 6:
         raise ValueError(f"Could not parse 6 header values from first line: '{header_str}'")

      # Metadata: ntot_u, ntot_v, ntot_c, T, dx, dy
      ntot_u, ntot_v, ntot_c = int(vec[0]), int(vec[1]), int(vec[2])
      T, dx, dy = vec[3], vec[4], vec[5]

      M_data = pd.read_csv(filepath, skiprows=1, header=None, usecols=range(7)).values

      return M_data, ntot_u, ntot_v, ntot_c, T, dx, dy

   except FileNotFoundError:
      print(f"  [ERROR] File not found: {filepath}")
      return None, 0, 0, 0, 0, 0, 0
   except Exception as e:
      print(f"  [ERROR] Processing {filepath}: {e}")
      return None, 0, 0, 0, 0, 0, 0


# --- 3. File setup arrays (0-indexed in Python) ---
ifile_s = 1
ifile_f = 3
ifile_ang = [0, np.arctan(1/2.), np.arctan(1.)]
ifile_str = ['0deg_', '27deg_', '45deg_']
jfile_str = ['stm', 'obs']

# List of dictionaries to replace MATLAB struct 'file'
files = []
skip_case = False

for ifile in range(ifile_s - 1, ifile_f): # indices 0, 1, 2
   jfile_s = 0
   jfile_f = 0
   if ifile == 0: # 0 deg case (stm and obs)
      jfile_f = 1
   for jfile in range(jfile_s, jfile_f + 1): # indices 0, 1 (for ifile 0), index 0 (for others)

      # Base file names for resolutions 32, 64, 128, 256, 320
      resolutions = [32, 64, 128, 256, 320]
      name_list = []
      for res in resolutions:
         name_list.append(
            f'rotated_cube_{ifile_str[ifile]}{res}_{jfile_str[jfile]}_mms.csv'
         )

      file_entry = {
         'name': name_list,
         'nameout': f'rotated_cube_{ifile_str[ifile]}{jfile_str[jfile]}_mms_convergence',
         'rotang': ifile_ang[ifile],
         'dx': [],
         'errors': {} # Dictionary to hold all computed error arrays
      }

      # Check existence and populate files list
      for name in file_entry['name']:
         if not os.path.exists(os.path.join(datadir, name)):
            skip_case = True
            print(f"[WARN] File {os.path.join(datadir, name)} does not exist. Skipping case.")

      files.append(file_entry)

if skip_case: quit()

# --- 4. Compute errors from *_mms.csv files ---
for ifile in range(len(files)):

   file_entry = files[ifile]
   rotangle = file_entry['rotang']

   ROTMAT = np.array([
      [np.cos(rotangle), -np.sin(rotangle)],
      [np.sin(rotangle), np.cos(rotangle)]
   ])
   TROTMAT = ROTMAT.T # Transpose

   # Initialize lists to store errors for all resolutions (n=1 to 5)
   error_lists = {k: [] for k in ['e_u_1', 'e_u_2', 'e_u_i', 'e_v_1', 'e_v_2', 'e_v_i',
                                 'e_z_1', 'e_z_2', 'e_z_i', 'e_H_1', 'e_H_2', 'e_H_i',
                                 'e_p_1', 'e_p_2', 'e_p_i']}

   for n in range(len(file_entry['name'])):
      filepath = os.path.join(datadir, file_entry['name'][n])

      M_data, ntot_u, ntot_v, ntot_c, T, dx, dy = read_matlab_csv(filepath)

      if M_data is None or M_data.size == 0:
         print(f"  Skipping processing for {file_entry['name'][n]} due to missing or empty data.")
         continue # Skip this resolution

      file_entry['dx'].append(dx)

      # Initialize error vectors
      e_u_vec = np.zeros(ntot_u)
      e_v_vec = np.zeros(ntot_v)
      e_z_vec = np.zeros(ntot_c)
      e_H_vec = np.zeros(ntot_c)
      e_p_vec = np.zeros(ntot_c)

      p = 0 # Data row counter (0-indexed)

      # --- First U: Velocity component along x-axis ---
      # M_data columns: 0: iscf, 1: xglob, 2: yglob, 3: area, 4: uglob
      for i in range(ntot_u):
         xglob = M_data[p, 1]
         yglob = M_data[p, 2]
         uglob = M_data[p, 4]

         # Compute analytical value in local coordinate system:
         XGLOB = np.array([xglob, yglob]) - displ
         XLOC = TROTMAT @ XGLOB - dispxy

         # Analytical Velocity in local axes (u, v):
         ULOC = np.array([
            ann_u(XLOC[0], XLOC[1], T),
            ann_v(XLOC[0], XLOC[1], T)
         ])

         # Analytical Velocity in global axes:
         UGLOB = ROTMAT @ ULOC
         uglob_ann = UGLOB[0]

         e_u_vec[i] = uglob - uglob_ann
         p += 1

      # Calculate norms for U
      error_lists['e_u_1'].append(np.linalg.norm(e_u_vec, 1) / ntot_u) # L1/N
      error_lists['e_u_2'].append(np.linalg.norm(e_u_vec, 2) / np.sqrt(ntot_u)) # L2/sqrt(N) (RMS-like)
      error_lists['e_u_i'].append(np.linalg.norm(e_u_vec, np.inf)) # L_inf


      # --- Second V: Velocity component along y-axis ---
      # M_data columns: 0: iscf, 1: xglob, 2: yglob, 3: area, 4: vglob
      for i in range(ntot_v):
         xglob = M_data[p, 1]
         yglob = M_data[p, 2]
         vglob = M_data[p, 4]

         # Compute analytical value in local coordinate system:
         XGLOB = np.array([xglob, yglob]) - displ
         XLOC = TROTMAT @ XGLOB - dispxy

         # Analytical Velocity in local axes (u, v):
         ULOC = np.array([
            ann_u(XLOC[0], XLOC[1], T),
            ann_v(XLOC[0], XLOC[1], T)
         ])

         # Analytical Velocity in global axes:
         UGLOB = ROTMAT @ ULOC
         vglob_ann = UGLOB[1]

         e_v_vec[i] = vglob - vglob_ann
         p += 1

      # Calculate norms for V
      error_lists['e_v_1'].append(np.linalg.norm(e_v_vec, 1) / ntot_v)
      error_lists['e_v_2'].append(np.linalg.norm(e_v_vec, 2) / np.sqrt(ntot_v))
      error_lists['e_v_i'].append(np.linalg.norm(e_v_vec, np.inf))

      # --- Finally cell centered variables (Z, H, P) ---
      # M_data columns: 0: iscf, 1: xglob, 2: yglob, 3: vol, 4: z, 5: H, 6: pres
      DELP = 0.0
      p_c = p # Start index for cell-centered variables

      for i in range(ntot_c):
         # Data from the CSV row:
         xglob = M_data[p, 1]
         yglob = M_data[p, 2]
         z = M_data[p, 4]
         H = M_data[p, 5]
         pres = M_data[p, 6]

         # Compute analytical scalars in local coordinate system:
         XGLOB = np.array([xglob, yglob]) - displ
         XLOC = TROTMAT @ XGLOB - dispxy

         z_ann = ann_z(XLOC[0], XLOC[1], T)
         H_ann = ann_H(XLOC[0], XLOC[1], T)
         p_ann = ann_p(XLOC[0], XLOC[1], T)

         if i == 0:
            # Pressure constant of integration (shift)
            DELP = pres - p_ann

         # Error calculation:
         e_z_vec[i] = z - z_ann
         e_H_vec[i] = H - H_ann - DELP / rho
         e_p_vec[i] = pres - p_ann - DELP # Pressure error is corrected by the shift DELP

         p += 1

      # Calculate norms for Z
      error_lists['e_z_1'].append(np.linalg.norm(e_z_vec, 1) / ntot_c)
      error_lists['e_z_2'].append(np.linalg.norm(e_z_vec, 2) / np.sqrt(ntot_c))
      error_lists['e_z_i'].append(np.linalg.norm(e_z_vec, np.inf))

      # Calculate norms for H
      error_lists['e_H_1'].append(np.linalg.norm(e_H_vec, 1) / ntot_c)
      error_lists['e_H_2'].append(np.linalg.norm(e_H_vec, 2) / np.sqrt(ntot_c))
      error_lists['e_H_i'].append(np.linalg.norm(e_H_vec, np.inf))

      # Calculate norms for P
      error_lists['e_p_1'].append(np.linalg.norm(e_p_vec, 1) / ntot_c)
      error_lists['e_p_2'].append(np.linalg.norm(e_p_vec, 2) / np.sqrt(ntot_c))
      error_lists['e_p_i'].append(np.linalg.norm(e_p_vec, np.inf))

   # Convert error lists to numpy arrays and store in the file entry
   for key, val_list in error_lists.items():
      file_entry['errors'][key] = np.array(val_list)

   files[ifile] = file_entry # Update the entry in the list


# --- 5. Generate Plots and Warnings ---
# print("\n--- Summary of Calculated Errors and Convergence Warnings ---")

for ifile in range(len(files)):
   file_entry = files[ifile]
   dx_values = np.array(file_entry['dx'])

   # Axis limits based on MATLAB original
   if ifile == 0 or ifile == 1:
      axisval=[1e-2, 1, 1e-7, 1e-1]
   elif ifile == 2:
      axisval=[1e-2, 1, 1e-6, 1e-1]
   else:
      axisval=[1e-2, 1, 1e-7, 1e-1]

   p_L2_u = np.log(file_entry['errors']['e_u_2'][:-1] / file_entry['errors']['e_u_2'][1:]) / np.log(dx_values[:-1] / dx_values[1:])
   # print(f"\nConvergence Order (p) for L2 Error u, {file_entry['nameout']}: {p_L2_u}")
   p_L2_z = np.log(file_entry['errors']['e_z_2'][:-1] / file_entry['errors']['e_z_2'][1:]) / np.log(dx_values[:-1] / dx_values[1:])
   # print(f"Convergence Order (p) for L2 Error z, {file_entry['nameout']}: {p_L2_z}")

   fig = fdsplotlib.plot_to_fig(x_data=dx_values, y_data=file_entry['errors']['e_z_2'], marker_style='rs-',plot_type='loglog',
         revision_label=version_string,x_min=axisval[0],x_max=axisval[1],y_min=axisval[2],y_max=axisval[3],
         data_label='FDS $Z$',
         x_label='$\\Delta x$ (m)',
         y_label='$L_2$ Error')

   fdsplotlib.plot_to_fig(x_data=dx_values, y_data= file_entry['errors']['e_u_2'], marker_style='g>-',
         figure_handle=fig,
         data_label='FDS $u$')

   fdsplotlib.plot_to_fig(x_data=dx_values, y_data= file_entry['errors']['e_H_2'], marker_style='b+>-',
         figure_handle=fig,
         data_label='FDS $H$')

   O1_ref = 10**-1 * dx_values
   fdsplotlib.plot_to_fig(x_data=dx_values, y_data= O1_ref, marker_style='k--',
         figure_handle=fig,
         data_label='$O(\\Delta x)$')
   O2_ref = 2 * 10**-3 * dx_values**2
   fdsplotlib.plot_to_fig(x_data=dx_values, y_data= O2_ref, marker_style='k-',
         figure_handle=fig,
         data_label='$O(\\Delta x^2)$')

   plotname = plotdir + file_entry['nameout'] + '.pdf'
   plt.savefig(plotname, format='pdf')
   plt.close()

   # --- Warning check logic ---
   last_z_error = file_entry['errors']['e_z_2'][-1]
   last_u_error = file_entry['errors']['e_u_2'][-1]
   last_H_error = file_entry['errors']['e_H_2'][-1]

   if ifile == 1: # OBST case (ifile == 2 in MATLAB)
      # print(f"\n*** Tolerance Check for OBST Case ({file_entry['nameout'].replace('_',' ')}) ***")
      if last_z_error > 2e-6:
         print(f"  Warning (Z): {last_z_error:.2e} > 2e-6 tolerance.")
      if last_u_error > 2e-5:
         print(f"  Warning (u): {last_u_error:.2e} > 2e-5 tolerance.")
      if last_H_error > 1.5e-3:
         print(f"  Warning (H): {last_H_error:.2e} > 1.5e-3 tolerance.")
   elif ifile == 2: # 45deg stm case (ifile == 3 in MATLAB)
      # print(f"\n*** Tolerance Check for 45deg STM Case ({file_entry['nameout'].replace('_',' ')}) ***")
      if last_z_error > 6.5e-5:
         print(f"  Warning (Z): {last_z_error:.2e} > 6.5e-5 tolerance.")
      if last_u_error > 8.5e-4:
         print(f"  Warning (u): {last_u_error:.2e} > 8.5e-4 tolerance.")
      if last_H_error > 2.5e-3:
         print(f"  Warning (H): {last_H_error:.2e} > 2.5e-3 tolerance.")


