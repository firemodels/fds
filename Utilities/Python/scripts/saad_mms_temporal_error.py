# McDermott
# 7-9-2015
# saad_mms_temporal_error.m
#
# Converted by Floyd
# 10-14-2025

import os
import numpy as np
import pandas as pd
import math
import matplotlib.pyplot as plt
import fdsplotlib


# Problem 1 parameters
r = 0.5
nx = 512
L = 2
dx = L/nx
x = np.arange(-L/2+dx/2, L/2, dx)
rho0 = 5
rho1 = .5
f = .5*(1+np.sin(2*np.pi*x/L))
rho = 1./((1-f)/rho0 + f/rho1)

# Output files

datadir = '../../Verification/Scalar_Analytical_Solution/';
plotdir = '../../Manuals/FDS_Verification_Guide/SCRIPT_FIGURES/';
filename = ['saad_512_cfl_1_mms.csv','saad_512_cfl_p5_mms.csv','saad_512_cfl_p25_mms.csv','saad_512_cfl_p125_mms.csv','saad_512_cfl_p0625_mms.csv']

skip_case = False

for i in range(len(filename)):
   name = datadir+filename[i]
   if not os.path.exists(name):
      skip_case = True
      print('Error: File ', filename[i], ' does not exist. Skipping case.')

if skip_case: quit()

# Gather FDS results

M1 = pd.read_csv(datadir+filename[2],skiprows=2,header=None)
M2 = pd.read_csv(datadir+filename[3],skiprows=2,header=None)
M3 = pd.read_csv(datadir+filename[4],skiprows=2,header=None)

ii = np.arange(nx,2*nx)
rho_1 = M1.iloc[ii,0]
rho_2 = M2.iloc[ii,0]
rho_3 = M3.iloc[ii,0]
Z_1 = M1.iloc[ii,1]
Z_2 = M2.iloc[ii,1]
Z_3 = M3.iloc[ii,1]

p_rho = np.log(np.abs(rho_3-rho_2)/np.abs(rho_2-rho_1) )/np.log(r);
p_Z = np.log( abs(Z_3-Z_2)/np.abs(Z_2-Z_1) )/np.log(r);

# print('Saad temporal order')
# print(' ')
# print('L1 p rho = ', np.linalg.norm(p_rho,ord=1)/nx)
# print('L2 p rho = ' ,np.linalg.norm(p_rho,2)/math.sqrt(nx))
# print('Linf p rho = ',np.linalg.norm(p_rho,ord=np.inf))
# print(' ')
# print('L1 p Z = ',np.linalg.norm(p_Z,ord=1)/nx)
# print('L2 p Z = ',np.linalg.norm(p_Z,ord=2)/math.sqrt(nx))
# print('Linf p Z = ',np.linalg.norm(p_Z,ord=npinf))
# print(' ')

# flag errors

L2_rho = np.linalg.norm(p_rho,ord=2)/math.sqrt(nx)
if L2_rho<1.99:
   print('Matlab Warning: L2_rho = ',L2_rho,' in Saad MMS')

L2_Z = np.linalg.norm(p_Z,ord=2)/math.sqrt(nx);
if L2_Z<1.99:
   print('Matlab Warning: L2_Z = ',L2_Z,' in Saad MMS')

# write the l2 norm to latex

with open(plotdir+'saad_l2_norm.tex', 'w') as fid:
   fid.write(f'{L2_rho:.4f}\n')

# Tony Saad's way...
p1_rho_saad = np.log( np.linalg.norm(rho_3-rho_2,ord=1)/np.linalg.norm(rho_2-rho_1,ord=1) )/np.log(r)
p2_rho_saad = np.log( np.linalg.norm(rho_3-rho_2,ord=2)/np.linalg.norm(rho_2-rho_1,ord=2) )/np.log(r)
pinf_rho_saad = np.log( np.linalg.norm(rho_3-rho_2,ord=np.inf)/np.linalg.norm(rho_2-rho_1,ord=np.inf) )/np.log(r)
# print('L1 p rho Saad = ', p1_rho_saad))
# print('L2 p rho Saad = ', p2_rho_saad )
# print('Linf p rho Saad = ', pinf_rho_saad))
# print(' ')
p1_Z_saad = np.log( np.linalg.norm(Z_3-Z_2,ord=1)/np.linalg.norm(Z_2-Z_1,ord=1) )/np.log(r)
p2_Z_saad = np.log( np.linalg.norm(Z_3-Z_2,ord=2)/np.linalg.norm(Z_2-Z_1,ord=2) )/np.log(r)
pinf_Z_saad = np.log( np.linalg.norm(Z_3-Z_2,ord=np.inf)/np.linalg.norm(Z_2-Z_1,ord=np.inf) )/np.log(r)
# print('L1 p Z Saad = ',p1_Z_saad )
# print('L2 p Z Saad = ',p2_Z_saad )
# print('Linf p Z Saad = ', pinf_Z_saad )
# print(' ')

git_file = datadir+'saad_512_cfl_p0625_git.txt'
version_string = fdsplotlib.get_version_string(git_file)

fig = fdsplotlib.plot_to_fig(x_data=x, y_data=p_rho, marker_style='b-',
      revision_label=version_string,x_min=-1,x_max=1,y_min=0,y_max=4,
      x_label='$x$ (m)',
      y_label='$p$ order density')

plotname = plotdir + 'saad_temporal_order_rho.pdf'
plt.savefig(plotname, format='pdf')
plt.close()

fig = fdsplotlib.plot_to_fig(x_data=x, y_data=p_Z, marker_style='b-',
      revision_label=version_string,x_min=-1,x_max=1,y_min=0,y_max=4,
      x_label='$x$ (m)',
      y_label='$p$ order mixture fraction')

plotname = plotdir + 'saad_temporal_order_Z.pdf'
plt.savefig(plotname, format='pdf')
plt.close()

fig = fdsplotlib.plot_to_fig(x_data=x, y_data=rho, marker_style='r--',
      revision_label=version_string,x_min=-1,x_max=1,y_min=0,y_max=5,
      data_label = 'Initial field',
      x_label='$x$ (m)',
      y_label='Density (kg/m$^3$')

fdsplotlib.plot_to_fig(x_data=x, y_data=rho_3, marker_style='r-',
   figure_handle=fig,
   data_label='Final field')

plotname = plotdir + 'saad_rho.pdf'
plt.savefig(plotname, format='pdf')
plt.close()

fig = fdsplotlib.plot_to_fig(x_data=x, y_data=f, marker_style='b--',
      revision_label=version_string,x_min=-1,x_max=1,y_min=0,y_max=1,
      data_label = 'Initial field',
      x_label='$x$ (m)',
      y_label='Mixture Fraction')

fdsplotlib.plot_to_fig(x_data=x, y_data=Z_3, marker_style='b-',
   figure_handle=fig,
   data_label='Final field')

plotname = plotdir + 'saad_Z.pdf'
plt.savefig(plotname, format='pdf')
plt.close()

fig = fdsplotlib.plot_to_fig(x_data=x, y_data=rho_3-rho_2, marker_style='b-',
      revision_label=version_string,x_min=-1,x_max=1,y_min=-0.0006,y_max=0.0006,
      data_label = '$\\rho_3$ - $\\rho_2$',
      x_label='$x$ (m)',
      y_label='Density (kg/m$^3$)')

fdsplotlib.plot_to_fig(x_data=x, y_data=rho_2-rho_1, marker_style='r-',
   figure_handle=fig,
      data_label = '$\\rho_2$ - $\\rho_1$')

plotname = plotdir + 'saad_rho_diff.pdf'
plt.savefig(plotname, format='pdf')
plt.close()

#Vanella - from McDermotts saad_mms_temporal_error.m
#02-20-2018
#saad_cc_mms_temporal_error.m
#
# Converted by Floyd
# 10-14-2025

datadir = '../../Verification/Complex_Geometry/'
filename = ['saad_CC_explicit_512_cfl_p25_mms.csv','saad_CC_explicit_512_cfl_p125_mms.csv','saad_CC_explicit_512_cfl_p0625_mms.csv']

skip_case = False

for i in range(len(filename)):
   name = datadir+filename[i]
   if not os.path.exists(name):
      skip_case = True
      print('Error: File ', filename[i], ' does not exist. Skipping case.')

if skip_case: quit()

M1 = pd.read_csv(datadir+filename[0],skiprows=2,header=None)
M2 = pd.read_csv(datadir+filename[1],skiprows=2,header=None)
M3 = pd.read_csv(datadir+filename[2],skiprows=2,header=None)

ii = np.arange(nx,2*nx)
rho_1 = M1.iloc[ii,0]
rho_2 = M2.iloc[ii,0]
rho_3 = M3.iloc[ii,0]
Z_1 = M1.iloc[ii,1]
Z_2 = M2.iloc[ii,1]
Z_3 = M3.iloc[ii,1]

p_rho = np.log(np.abs(rho_3-rho_2)/np.abs(rho_2-rho_1) )/np.log(r);
p_Z = np.log( abs(Z_3-Z_2)/np.abs(Z_2-Z_1) )/np.log(r);

# print('Saad CC temporal order')
# print(' ')
# print('L1 p rho = ', np.linalg.norm(p_rho,ord=1)/nx)
# print('L2 p rho = ' ,np.linalg.norm(p_rho,2)/math.sqrt(nx))
# print('Linf p rho = ',np.linalg.norm(p_rho,ord=np.inf))
# print(' ')
# print('L1 p Z = ',np.linalg.norm(p_Z,ord=1)/nx)
# print('L2 p Z = ',np.linalg.norm(p_Z,ord=2)/math.sqrt(nx))
# print('Linf p Z = ',np.linalg.norm(p_Z,ord=npinf))
# print(' ')

# flag errors

L2_rho = np.linalg.norm(p_rho,ord=2)/math.sqrt(nx)
if L2_rho<1.99:
   print('Matlab Warning: L2_rho = ',L2_rho,' in Saad CC MMS')

L2_Z = np.linalg.norm(p_Z,ord=2)/math.sqrt(nx);
if L2_Z<1.99:
   print('Matlab Warning: L2_Z = ',L2_Z,' in Saad CC MMS')

# write the l2 norm to latex

with open(plotdir+'saad_CC_explicit_l2_norm.tex', 'w') as fid:
   fid.write(f'{L2_rho:.4f}\n')

# Tony Saad's way...
p1_rho_saad = np.log( np.linalg.norm(rho_3-rho_2,ord=1)/np.linalg.norm(rho_2-rho_1,ord=1) )/np.log(r)
p2_rho_saad = np.log( np.linalg.norm(rho_3-rho_2,ord=2)/np.linalg.norm(rho_2-rho_1,ord=2) )/np.log(r)
pinf_rho_saad = np.log( np.linalg.norm(rho_3-rho_2,ord=np.inf)/np.linalg.norm(rho_2-rho_1,ord=np.inf) )/np.log(r)
# print('L1 p rho Saad = ', p1_rho_saad))
# print('L2 p rho Saad = ', p2_rho_saad )
# print('Linf p rho Saad = ', pinf_rho_saad))
# print(' ')
p1_Z_saad = np.log( np.linalg.norm(Z_3-Z_2,ord=1)/np.linalg.norm(Z_2-Z_1,ord=1) )/np.log(r)
p2_Z_saad = np.log( np.linalg.norm(Z_3-Z_2,ord=2)/np.linalg.norm(Z_2-Z_1,ord=2) )/np.log(r)
pinf_Z_saad = np.log( np.linalg.norm(Z_3-Z_2,ord=np.inf)/np.linalg.norm(Z_2-Z_1,ord=np.inf) )/np.log(r)
# print('L1 p Z Saad = ',p1_Z_saad )
# print('L2 p Z Saad = ',p2_Z_saad )
# print('Linf p Z Saad = ', pinf_Z_saad )
# print(' ')

git_file = datadir+'saad_CC_explicit_512_cfl_p0625_git.txt'
version_string = fdsplotlib.get_version_string(git_file)

fig = fdsplotlib.plot_to_fig(x_data=x, y_data=p_rho, marker_style='b-',
      revision_label=version_string,x_min=-1,x_max=1,y_min=0,y_max=4,
      x_label='$x$ (m)',
      y_label='$p$ order density')

plotname = plotdir + 'saad_CC_explicit_temporal_order_rho.pdf'
plt.savefig(plotname, format='pdf')
plt.close()

fig = fdsplotlib.plot_to_fig(x_data=x, y_data=p_Z, marker_style='b-',
      revision_label=version_string,x_min=-1,x_max=1,y_min=0,y_max=4,
      x_label='$x$ (m)',
      y_label='$p$ order mixture fraction')

plotname = plotdir + 'saad_CC_explicit_temporal_order_Z.pdf'
plt.savefig(plotname, format='pdf')
plt.close()

fig = fdsplotlib.plot_to_fig(x_data=x, y_data=rho, marker_style='r--',
      revision_label=version_string,x_min=-1,x_max=1,y_min=0,y_max=5,
      data_label = 'Initial field',
      x_label='$x$ (m)',
      y_label='Density (kg/m$^3$')

fdsplotlib.plot_to_fig(x_data=x, y_data=rho_3, marker_style='r-',
   figure_handle=fig,
   data_label='Final field')

plotname = plotdir + 'saad_CC_explicit_rho.pdf'
plt.savefig(plotname, format='pdf')
plt.close()

fig = fdsplotlib.plot_to_fig(x_data=x, y_data=f, marker_style='b--',
      revision_label=version_string,x_min=-1,x_max=1,y_min=0,y_max=1,
      data_label = 'Initial field',
      x_label='$x$ (m)',
      y_label='Mixture Fraction')

fdsplotlib.plot_to_fig(x_data=x, y_data=Z_3, marker_style='b-',
   figure_handle=fig,
   data_label='Final field')

plotname = plotdir + 'saad_CC_explicit_Z.pdf'
plt.savefig(plotname, format='pdf')
plt.close()

fig = fdsplotlib.plot_to_fig(x_data=x, y_data=rho_3-rho_2, marker_style='b-',
      revision_label=version_string,x_min=-1,x_max=1,y_min=-0.0006,y_max=0.0006,
      data_label = '$\\rho_3$ - $\\rho_2$',
      x_label='$x$ (m)',
      y_label='Density (kg/m$^3$)')

fdsplotlib.plot_to_fig(x_data=x, y_data=rho_2-rho_1, marker_style='r-',
   figure_handle=fig,
      data_label = '$\\rho_2$ - $\\rho_1$')

plotname = plotdir + 'saad_CC_explicit_rho_diff.pdf'
plt.savefig(plotname, format='pdf')
plt.close()
