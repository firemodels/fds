# Roger Wang
# 6-28-11
# ns2d.py
# Process the output of the ns2d verification cases
#
# Converted by Floyd
# 10/17/2025

import numpy as np
import pandas as pd
import os
import math


datadir = '../../Verification/NS_Analytical_Solution/'
infile = [['ns2d_8_devc.csv','ns2d_8_nupt1_devc.csv'],
   ['ns2d_16_devc.csv','ns2d_16_nupt1_devc.csv'],
   ['ns2d_32_devc.csv','ns2d_32_nupt1_devc.csv'],
   ['ns2d_64_devc.csv','ns2d_64_nupt1_devc.csv']]

outfile = [['ns2d_8_exact.csv','ns2d_8_nupt1_exact.csv'],
   ['ns2d_16_exact.csv','ns2d_16_nupt1_exact.csv'],
   ['ns2d_32_exact.csv','ns2d_32_nupt1_exact.csv'],
   ['ns2d_64_exact.csv','ns2d_64_nupt1_exact.csv']]

skip_case = False

for i in range(len(infile)):
   for j in range(len(infile[i])):
      name = datadir+infile[i][j]
      if not os.path.exists(name):
         skip_case = True
         print('Error: File ', infile[i][j], ' does not exist. Skipping case.')

if skip_case: quit()

nu = [0.0, 0.1]
pi = math.pi
x = pi
y = pi
dx = [2.0 * pi / 8.0, 2.0 * pi / 16.0, 2.0 * pi / 32.0, 2.0 * pi / 64.0]

rms = np.zeros((4, 2))
cnt = np.zeros((4, 2))

# k iterates over viscosity (0: nu=0.0, 1: nu=0.1)
for k in range(2):
   # j iterates over resolution (0: 8 cells, 3: 64 cells)
   for j in range(4):
      M_10 = pd.read_csv(datadir+infile[j][k], skiprows=2, header=None).values

      with open(datadir+outfile[j][k], 'w') as M_11:
         M_11.write("Time, u-vel\n")

         squared_error_sum = 0.0
         num_points = 0.0

         for i in range(len(M_10)):
            t = M_10[i, 0]
            u_fds = M_10[i, 1]

            u_exact = 1.0 - 2.0 * math.cos(x - t) * math.sin(y - 0.5 * dx[j] - t) * math.exp(-2.0 * nu[k] * t)

            M_11.write(f"{t:8.3f}, {u_exact:9.5f}\n")

            squared_error_sum += (u_fds - u_exact)**2

            num_points += 1.0

         if num_points > 0:
            rms[j, k] = math.sqrt(squared_error_sum / num_points)
         else:
            rms[j, k] = 0.0

# Case 1: nu=0.0 (k=0)
with open(datadir+'ns2d_error.csv', 'w') as fid11:
   fid11.write('dx,dx^2,rms error\n')
   for j in range(4):
      fid11.write(f"{dx[j]:8.3f}, {dx[j]**2:8.3f}, {rms[j, 0]:8.4f}\n")

# Case 2: nu=0.1 (k=1)
with open(datadir+'ns2d_nupt1_error.csv', 'w') as fid12:
   fid12.write('dx,dx^2,rms error\n')
   for j in range(4):
      fid12.write(f"{dx[j]:8.3f}, {dx[j]**2:8.3f}, {rms[j, 1]:8.4f}\n")

