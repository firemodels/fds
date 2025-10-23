# shun_mms_py
#!/usr/bin/env python3
# McDermott
# 8-3-2015
# shunn_mms_temporal_error.m
#
# Converted by Floyd
# 10-14-2025

import os
import numpy as np
import pandas as pd
import math
import matplotlib.pyplot as plt
import fdsplotlib

r =0.5
nx = 256

datadir = '../../Verification/Scalar_Analytical_Solution/'
plotdir = '../../Manuals/FDS_Verification_Guide/SCRIPT_FIGURES/'
filename = ['shunn3_256_cfl_1_mms.csv',
   'shunn3_256_cfl_p5_mms.csv',
   'shunn3_256_cfl_p25_mms.csv',
   'shunn3_256_cfl_p125_mms.csv',
   'shunn3_256_cfl_p0625_mms.csv']

skip_case = False

for i in range(len(filename)):
   name = datadir+filename[i]
   if not os.path.exists(name):
      skip_case = True
      print('Error: File ', filename[i], ' does not exist. Skipping case.')

if skip_case: quit()

M1 = pd.read_csv(datadir+filename[2],skiprows=2,header=None)
M2 = pd.read_csv(datadir+filename[3],skiprows=2,header=None)
M3 = pd.read_csv(datadir+filename[4],skiprows=2,header=None)

rho_1 = M1.iloc[:,0]
rho_2 = M2.iloc[:,0]
rho_3 = M3.iloc[:,0]
Z_1 = M1.iloc[:,1]
Z_2 = M2.iloc[:,1]
Z_3 = M3.iloc[:,1]
U_1 = M1.iloc[:,2]
U_2 = M2.iloc[:,2]
U_3 = M3.iloc[:,2]
H_1 = M1.iloc[:,4]
H_2 = M2.iloc[:,4]
H_3 = M3.iloc[:,4]

p_rho = np.log(np.abs(rho_3-rho_2)/np.abs(rho_2-rho_1) )/np.log(r);
p_Z = np.log(np.abs(Z_3-Z_2)/np.abs(Z_2-Z_1) )/np.log(r);
p_U = np.log(np.abs(U_3-U_2)/np.abs(U_2-U_1) )/np.log(r);
p_H = np.log(np.abs(H_3-H_2)/np.abs(H_2-H_1) )/np.log(r);


# print('Shunn 3 temporal order')
# print(' ')
# print('L1 p rho = ',( np.linalg.norm(p_rho,ord=1)/(nx*nx) ))
# print('L2 p rho = ',( np.linalg.norm(p_rho,2ord=)/nx ))
# print('Linf p rho = ',( np.linalg.norm(p_rho,ord=np.inf) ))
# print(' ')
# print('L1 p Z = ',( np.linalg.norm(p_Z,ord=1)/(nx*nx) ))
# print('L2 p Z = ',( np.linalg.norm(p_Z,ord=2)/nx ))
# print('Linf p Z = ',( np.linalg.norm(p_Z,ord=np.inf) ))
# print(' ')
# print('L1 p U = ',( np.linalg.norm(p_U,ord=1)/(nx*nx) ))
# print('L2 p U = ',( np.linalg.norm(p_U,ord=2)/nx ))
# print('Linf p U = ',( np.linalg.norm(p_U,ord=np.nf) ))
# print(' ')
# print('L1 p H = ',( np.linalg.norm(p_H,ord=1)/(nx*nx) ))
# print('L2 p H = ',( np.linalg.norm(p_H,ord=2)/nx ))
# print('Linf p H = ',( np.linalg.norm(p_H,ord=np.inf) ))
# print(' ')

L2_rho = np.linalg.norm(p_rho,ord=2)/nx;
if L2_rho<1.99:
    disp(['Python Warning: L2_rho = ',L2_rho,' in Shunn 3 MMS temporal order'])

L2_Z = np.linalg.norm(p_Z,ord=2)/nx;
if L2_Z<1.99:
   print('Python Warning: L2_Z = ',L2_Z,' in Shunn 3 MMS temporal order')

L2_U = np.linalg.norm(p_U,ord=2)/nx;
if L2_U<1.99:
   print('Python Warning: L2_U = ',L2_U,' in Shunn 3 MMS temporal order')

L2_H = np.linalg.norm(p_H,ord=2)/nx;
if L2_H<0.99:
   print('Python Warning: L2_H = ',L2_H,' in Shunn 3 MMS temporal order')


# Tony Saad's way...
p1_rho_saad = np.log( np.linalg.norm(rho_3-rho_2,1)/np.linalg.norm(rho_2-rho_1,ord=1) )/np.log(r);
p2_rho_saad = np.log( np.linalg.norm(rho_3-rho_2,2)/np.linalg.norm(rho_2-rho_1,ord=2) )/np.log(r);
pinf_rho_saad = np.log( np.linalg.norm(rho_3-rho_2,ord=np.inf)/np.linalg.norm(rho_2-rho_1,ord=np.inf) )/np.log(r);
# print('L1 p rho Saad = ',num2str( p1_rho_saad ))
# print('L2 p rho Saad = ',num2str( p2_rho_saad ))
# print('Lord=np.inf p rho Saad = ',num2str( pord=np.inf_rho_saad ))
# print(' ')
p1_Z_saad = np.log( np.linalg.norm(Z_3-Z_2,1)/np.linalg.norm(Z_2-Z_1,ord=1) )/np.log(r);
p2_Z_saad = np.log( np.linalg.norm(Z_3-Z_2,2)/np.linalg.norm(Z_2-Z_1,ord=2) )/np.log(r);
pinf_Z_saad = np.log( np.linalg.norm(Z_3-Z_2,ord=np.inf)/np.linalg.norm(Z_2-Z_1,ord=np.inf) )/np.log(r);
# print('L1 p Z Saad = ',num2str( p1_Z_saad ))
# print('L2 p Z Saad = ',num2str( p2_Z_saad ))
# print('Lord=np.inf p Z Saad = ',num2str( pord=np.inf_Z_saad ))
# print(' ')
p1_U_saad = np.log( np.linalg.norm(U_3-U_2,1)/np.linalg.norm(U_2-U_1,ord=1) )/np.log(r);
p2_U_saad = np.log( np.linalg.norm(U_3-U_2,2)/np.linalg.norm(U_2-U_1,ord=2) )/np.log(r);
pinf_U_saad = np.log( np.linalg.norm(U_3-U_2,ord=np.inf)/np.linalg.norm(U_2-U_1,ord=np.inf) )/np.log(r);
# print('L1 p U Saad = ',num2str( p1_U_saad ))
# print('L2 p U Saad = ',num2str( p2_U_saad ))
# print('Lord=np.inf p U Saad = ',num2str( pord=np.inf_U_saad ))
# print(' ')
p1_H_saad = np.log( np.linalg.norm(H_3-H_2,1)/np.linalg.norm(H_2-H_1,ord=1) )/np.log(r);
p2_H_saad = np.log( np.linalg.norm(H_3-H_2,2)/np.linalg.norm(H_2-H_1,ord=2) )/np.log(r);
pinf_H_saad = np.log( np.linalg.norm(H_3-H_2,ord=np.inf)/np.linalg.norm(H_2-H_1,ord=np.inf) )/np.log(r);
# print('L1 p H Saad = ',num2str( p1_H_saad ))
# print('L2 p H Saad = ',num2str( p2_H_saad ))
# print('Linf p H Saad = ',num2str( pinf_H_saad ))

# write the norms to a latex table

with open(plotdir+'shunn_terr_norms.tex', 'w') as fid:
   tmp_str = f'$L_1$ & {np.linalg.norm(p_rho, ord=1) / (nx * nx):.4f} & {np.linalg.norm(p_Z, ord=1) / (nx * nx):.4f} & {np.linalg.norm(p_U, ord=1) / (nx * nx):.4f} & {np.linalg.norm(p_H, ord=1) / (nx * nx):.4f} \\\\'
   fid.write(f"{tmp_str}\n")
   tmp_str = f'$L_2$ & {np.linalg.norm(p_rho, ord=2) / (nx * nx):.4f} & {np.linalg.norm(p_Z, ord=2) / (nx * nx):.4f} & {np.linalg.norm(p_U, ord=2) / (nx * nx):.4f} & {np.linalg.norm(p_H, ord=2) / (nx * nx):.4f} \\\\'
   fid.write(f"{tmp_str}\n")
   tmp_str = f'$L_1$ Saad & {p1_rho_saad:.4f} & {p1_Z_saad:.4f} & {p1_U_saad:.4f} & {p1_H_saad:.4f} \\\\'
   fid.write(f"{tmp_str}\n")
   tmp_str = f'$L_2$ Saad & {p2_rho_saad:.4f} & {p2_Z_saad:.4f} & {p2_U_saad:.4f} & {p2_H_saad:.4f} \\\\'
   fid.write(f"{tmp_str}\n")
   tmp_str = f'$L_{{\\infty}}$ Saad & {pinf_rho_saad:.4f} & {pinf_Z_saad:.4f} & {pinf_U_saad:.4f} & {pinf_H_saad:.4f} \\\\'
   fid.write(f"{tmp_str}\n")

# McDermott
# 9-4-2013
# shunn_mms_favreZ.m
#
# Converted by Floyd
# 10-14-2025

r0 = 5.
r1 = 1.
uf = 0.5
vf = 0.5
k = 2.
w = 2.
mu = 0.001
D  = 0.001


def vd2d_mms_z(x, y, t):
   numerator = 1.0 + math.sin(math.pi * k * (x - uf * t)) * \
   math.sin(math.pi * k * (y - vf * t)) * \
   math.cos(math.pi * w * t)

   denominator = (1 + r0 / r1) + (1 - r0 / r1) * \
   math.sin(math.pi * k * (x - uf * t)) * \
   math.sin(math.pi * k * (y - vf * t)) * \
   math.cos(math.pi * w * t)

   return numerator / denominator

def vd2d_mms_rho(x, y, t):
   numerator = 1.0
   denominator =vd2d_mms_z(x, y, t)/r1 + (1-vd2d_mms_z(x, y, t))/r0

   return numerator / denominator

X = 0
Y = 0
T = 10
N = 1000
dt = T/N

Zbar    = 0
rhoZbar = 0
rhobar  = 0
for n in list(range(1,N+1)):
   t = n*dt
   Zbar    = Zbar    + vd2d_mms_z(X,Y,t)*dt
   rhoZbar = rhoZbar + vd2d_mms_rho(X,Y,t)*vd2d_mms_z(X,Y,t)*dt
   rhobar  = rhobar  + vd2d_mms_rho(X,Y,t)*dt

Zbar_mms = Zbar/T
rhoZbar = rhoZbar/T
rhobar = rhobar/T

FavreZ_mms = rhoZbar/rhobar

filename = ['shunn3_FavreZ_32_devc.csv','shunn3_FavreZ_64_devc.csv']
git_file = datadir+'shunn3_FavreZ_32_git.txt'
version_string = fdsplotlib.get_version_string(git_file)

skip_case = False

for i in range(len(filename)):
   name = datadir+filename[i]
   if not os.path.exists(name):
      skip_case = True
      print('Error: File ', filename[i], ' does not exist. Skipping case.')

if skip_case: quit()

fig = fdsplotlib.plot_to_fig(x_data=[0,T], y_data=[FavreZ_mms,FavreZ_mms], marker_style='k-',
      revision_label=version_string,x_min=0,x_max=T,y_min=0,y_max=0.2,
      data_label='Analytical Favre',
      x_label='Time (s)',
      y_label='Mass Fraction')

fdsplotlib.plot_to_fig(x_data=[0,T], y_data=[Zbar_mms,Zbar_mms], marker_style='k--',
      figure_handle=fig,
      data_label='Analytical Average')

linestyle = ['m-','m--','b-','b--']
pltlabel = ['FDS 32 Favre','FDS 32 Average','FDS 64 Favre','FDS 64 Average']

k=0
for i in range(len(filename)):
   M = pd.read_csv(datadir+filename[i],skiprows=2,header=None)
   t = M.iloc[:,0]
   Zbar = M.iloc[:,1]
   FavreZ= M.iloc[:,2]
   fdsplotlib.plot_to_fig(x_data=t, y_data=FavreZ, marker_style=linestyle[k],
      figure_handle=fig,
      data_label=pltlabel[k])
   k = k + 1
   fdsplotlib.plot_to_fig(x_data=t, y_data=Zbar, marker_style=linestyle[k],
      figure_handle=fig,
      data_label=pltlabel[k])
   k = k + 1

plotname = plotdir + 'shunn_mms_FavreZ.pdf'
plt.savefig(plotname, format='pdf')
plt.close()

# check errors

FavreZ_error_64 = np.abs(FavreZ[len(FavreZ)-1]-FavreZ_mms)/FavreZ_mms;
if FavreZ_error_64 > 0.01:
   print('Python Warning: FavreZ in shunn3_FavreZ_64 is out of tolerance. Error = ',FavreZ_error_64)


#McDermott
#09-04-2013
# shunn_mms_error.m
#
# Converted by Floyd
# 10-15-2025

r0 = 5.
r1 = 1.
uf = 0.5
vf = 0.5
k = 2.
w = 2.
mu = 0.001
D  = 0.001

# analytical solutions
def vd2d_mms_z (x,y,t) :
   numerator = ( 1. + math.sin(math.pi*k*(x-uf*t))*math.sin(math.pi*k*(y-vf*t))*math.cos(math.pi*w*t) )
   denominator = ( (1+r0/r1) + (1-r0/r1)*math.sin(math.pi*k*(x-uf*t))*math.sin(math.pi*k*(y-vf*t))*math.cos(math.pi*w*t) )
   return numerator/denominator

def vd2d_mms_rho(x,y,t):
   value=1./( vd2d_mms_z(x,y,t)/r1 + (1-vd2d_mms_z(x,y,t))/r0 )
   return value

def vd2d_mms_u(x,y,t):
   value = uf + (r1-r0)/vd2d_mms_rho(x,y,t)*(-w/(4.*k))*math.cos(math.pi*k*(x-uf*t))*math.sin(math.pi*k*(y-vf*t))*math.sin(math.pi*w*t)
   return value

def vd2d_mms_v(x,y,t):
   value = vf + (r1-r0)/vd2d_mms_rho(x,y,t)*(-w/(4.*k))*math.sin(math.pi*k*(x-uf*t))*math.cos(math.pi*k*(y-vf*t))*math.sin(math.pi*w*t)
   return value

def vd2d_mms_H(x,y,t):
   value = 0.5*(vd2d_mms_u(x,y,t)-uf)*(vd2d_mms_v(x,y,t)-vf)
   return value

def vd2d_mms_p(x,y,t) :
   vd2d_mms_rho(x,y,t)*(vd2d_mms_H(x,y,t) - 0.5*(vd2d_mms_u(x,y,t)**2 + vd2d_mms_v(x,y,t)**2))
   return value

L = 2
nx = np.array([32,64,128,256,512])
dx = L/nx

datadir = '../../Verification/Scalar_Analytical_Solution/';
filename = ['shunn3_32_mms.csv','shunn3_64_mms.csv','shunn3_128_mms.csv','shunn3_256_mms.csv','shunn3_512_mms.csv']

skip_case = False

for i in range(len(filename)):
   name = datadir+filename[i]
   if not os.path.exists(name):
      skip_case = True
      print('Error: File ', filename[i], ' does not exist. Skipping case.')

if skip_case: quit()

e_r = np.zeros(len(filename))
e_z = np.zeros(len(filename))
e_u = np.zeros(len(filename))
e_H = np.zeros(len(filename))

for n in range(len(filename)):
   N = nx[n]**2


   #print(filename[n])
   M = pd.read_csv(datadir+filename[n],skiprows=1,header=None,names=['1','2','3','4','5','6'])

   T = M.iloc[0,0]

   M = pd.read_csv(datadir+filename[n],skiprows=2,header=None)

   x = np.linspace(-L/2, L/2, nx[n] + 1)
   xc = x[:nx[n]]+0.5*dx[n]
   yc = xc

   #inialize error arrays
   rho_error = np.zeros((nx[n],nx[n]))

   e_r_vec = np.zeros(N)

   z_error = np.zeros((nx[n],nx[n]))
   e_z_vec = np.zeros(N)

   u_error = np.zeros((nx[n],nx[n]))
   e_u_vec = np.zeros(N)

   H_error = np.zeros((nx[n],nx[n]))
   e_H_vec = np.zeros(N)

   p = 0;
   for j in range(nx[n]):
      for i in range(nx[n]):
         rho = M.iloc[p,0]
         rho_mms = vd2d_mms_rho(xc[i],yc[j],T)
         rho_error[i,j] = rho - rho_mms
         e_r_vec[p] = rho_error[i,j]

         z = M.iloc[p,1]
         z_mms = vd2d_mms_z(xc[i],yc[j],T)
         z_error[i,j] = z - z_mms
         e_z_vec[p] = z_error[i,j]

         u = M.iloc[p,2]
         u_mms = vd2d_mms_u(x[i+1],yc[j],T)
         u_error[i,j] = u - u_mms
         e_u_vec[p] = u_error[i,j]

         H = M.iloc[p,4]
         H_mms = vd2d_mms_H(xc[i],yc[j],T)
         H_error[i,j] = H - H_mms
         e_H_vec[p] = H_error[i,j]

         p = p+1

   e_r[n] = np.linalg.norm(e_r_vec,ord=2)/nx[n]
   e_z[n] = np.linalg.norm(e_z_vec,ord=2)/nx[n]
   e_u[n] = np.linalg.norm(e_u_vec,ord=2)/nx[n]
   e_H[n] = np.linalg.norm(e_H_vec,ord=2)/nx[n]

# add Git version if file is available
git_file = datadir+'shunn3_256_git.txt'
version_string = fdsplotlib.get_version_string(git_file)

datalabel=['FDS $\\rho$','FDS $Z$','FDS $u$','FDS $H$','$O$($\\Delta x$)','$O$($\\Delta x^2$)']

fig = fdsplotlib.plot_to_fig(x_data=dx, y_data=e_r, marker_style='ro-',plot_type='loglog',
      revision_label=version_string,x_min=0.001,x_max=0.1,y_min=1E-6,y_max=0.1,
      data_label=datalabel[0],
      x_label='$\\Delta x$ (m)',
      y_label='$L_2$ Error')

fdsplotlib.plot_to_fig(x_data=dx, y_data=e_z, marker_style='gs-',
      figure_handle=fig,
      data_label=datalabel[1])

fdsplotlib.plot_to_fig(x_data=dx, y_data=e_u, marker_style='b>-',
      figure_handle=fig,
      data_label=datalabel[2])

fdsplotlib.plot_to_fig(x_data=dx, y_data=e_H, marker_style='c+-',
      figure_handle=fig,
      data_label=datalabel[3])

fdsplotlib.plot_to_fig(x_data=dx, y_data=dx, marker_style='k--',
      figure_handle=fig,
      data_label=datalabel[4])

fdsplotlib.plot_to_fig(x_data=dx, y_data=dx**2, marker_style='k-',
      figure_handle=fig,
      data_label=datalabel[5])

plotname = plotdir + 'shunn_mms_convergence.pdf'
plt.savefig(plotname, format='pdf')
plt.close()

# check errors
if e_r[len(e_r)-1] > 2e-4:
   print('Python Warning: Density in shunn3 is out of tolerance. e_r = ',e_r[len(e_r)-1] )
if e_z[len(e_z)-1] > 1e-4:
   print('Python Warning: Mixture fraction in shunn3 is out of tolerance. e_z = ',e_z[len(e_Z)-1] )
if e_u[len(e_u)-1] > 3e-5:
   print('Python Warning: Velocity in shunn3 is out of tolerance. e_u = ',e_u[len(e_u)-1] )
if e_H[len(e_H)-1] > 2e-3:
   print('Python Warning: Pressure in shunn3 is out of tolerance. e_H = ',e_H[len(e_H)-1] )

#Vanella from McDermott
#03-09-2018
# shunn_CC_mms_error.m
#
# Converted by Floyd
# 10-15-2025

nx = np.array([32,64,128,256,320])
dx = L/nx
ymax=[0.1,1]

datadir = '../../Verification/Complex_Geometry/'
filename = [['shunn3_32_cc_exp_chm_mms.csv','shunn3_64_cc_exp_chm_mms.csv','shunn3_128_cc_exp_chm_mms.csv',
   'shunn3_256_cc_exp_chm_mms.csv','shunn3_320_cc_exp_chm_mms.csv'],
   ['shunn3_32_cc_exp_gdv_mms.csv','shunn3_64_cc_exp_gdv_mms.csv','shunn3_128_cc_exp_gdv_mms.csv','shunn3_256_cc_exp_gdv_mms.csv',
   'shunn3_320_cc_exp_gdv_mms.csv']]
gitname = ['shunn3_256_cc_exp_chm_git.txt','shunn3_256_cc_exp_gdv_git.txt']
plotname=[ 'shunn_cc_exp_chm_mms_convergence.pdf', 'shunn_cc_exp_gdv_mms_convergence.pdf']

skip_case = False

for i in range(2):
   for j in range(len(filename[i])):
      name = datadir+filename[i][j]
      if not os.path.exists(name):
         skip_case = True
         print('Error: File ', filename[i][j], ' does not exist. Skipping case.')

if skip_case: quit()

for kk in range(2):
   for n in range(len(filename[kk])):
      N = nx[n]**2

      #print(filename[n])
      M = pd.read_csv(datadir+filename[kk][n],skiprows=1,header=None,names=['1','2','3','4','5','6'])

      T = M.iloc[0,0]

      M = pd.read_csv(datadir+filename[kk][n],skiprows=2,header=None)

      x = np.linspace(-L/2, L/2, nx[n] + 1)
      xc = x[:nx[n]]+0.5*dx[n]
      yc = xc

      #inialize error arrays
      rho_error = np.zeros((nx[n],nx[n]))

      e_r_vec = np.zeros(N)

      z_error = np.zeros((nx[n],nx[n]))
      e_z_vec = np.zeros(N)

      u_error = np.zeros((nx[n],nx[n]))
      e_u_vec = np.zeros(N)

      H_error = np.zeros((nx[n],nx[n]))
      e_H_vec = np.zeros(N)

      p = 0;
      for j in range(nx[n]):
         for i in range(nx[n]):
            rho = M.iloc[p,0]
            rho_mms = vd2d_mms_rho(xc[i],yc[j],T)
            rho_error[i,j] = rho - rho_mms
            e_r_vec[p] = rho_error[i,j]

            z = M.iloc[p,1]
            z_mms = vd2d_mms_z(xc[i],yc[j],T)
            z_error[i,j] = z - z_mms
            e_z_vec[p] = z_error[i,j]

            u = M.iloc[p,2]
            u_mms = vd2d_mms_u(x[i+1],yc[j],T)
            u_error[i,j] = u - u_mms
            e_u_vec[p] = u_error[i,j]

            H = M.iloc[p,4]
            H_mms = vd2d_mms_H(xc[i],yc[j],T)
            H_error[i,j] = H - H_mms
            e_H_vec[p] = H_error[i,j]

            p = p+1

      e_r[n] = np.linalg.norm(e_r_vec,ord=2)/nx[n]
      e_z[n] = np.linalg.norm(e_z_vec,ord=2)/nx[n]
      e_u[n] = np.linalg.norm(e_u_vec,ord=2)/nx[n]
      e_H[n] = np.linalg.norm(e_H_vec,ord=2)/nx[n]

      # add Git version if file is available
   git_file = datadir+gitname[kk]
   version_string = fdsplotlib.get_version_string(git_file)

   datalabel=['FDS $\\rho$','FDS $Z$','FDS $u$','FDS $H$','$O$($\\Delta x$)','$O$($\\Delta x^2$)']

   fig = fdsplotlib.plot_to_fig(x_data=dx, y_data=e_r, marker_style='ro-',plot_type='loglog',
         revision_label=version_string,x_min=0.001,x_max=0.1,y_min=1E-5,y_max=ymax[kk],
         data_label=datalabel[0],
         x_label='$\\Delta x$ (m)',
         y_label='$L_2$ Error')

   fdsplotlib.plot_to_fig(x_data=dx, y_data=e_z, marker_style='gs-',
         figure_handle=fig,
         data_label=datalabel[1])

   fdsplotlib.plot_to_fig(x_data=dx, y_data=e_u, marker_style='b>-',
         figure_handle=fig,
         data_label=datalabel[2])

   fdsplotlib.plot_to_fig(x_data=dx, y_data=e_H, marker_style='c+-',
         figure_handle=fig,
         data_label=datalabel[3])

   fdsplotlib.plot_to_fig(x_data=dx, y_data=dx, marker_style='k--',
         figure_handle=fig,
         data_label=datalabel[4])

   fdsplotlib.plot_to_fig(x_data=dx, y_data=dx**2, marker_style='k-',
         figure_handle=fig,
         data_label=datalabel[5])

   plotfile = plotdir + plotname[kk]
   plt.savefig(plotfile, format='pdf')
   plt.close()

   if(kk==0):
      # check errors
      if e_r[len(e_r)-1] > 3e-4:
         print('Python Warning: Density in shunn3 CC is out of tolerance. e_r = ',e_r[len(e_r)-1] )
      if e_z[len(e_z)-1] > 1e-4:
         print('Python Warning: Mixture fraction in shunn3 CC is out of tolerance. e_z = ',e_z[len(e_z)-1] )
      if e_u[len(e_u)-1] > 3e-5:
         print('Python Warning: Velocity in shunn3 CC is out of tolerance. e_u = ',e_u[len(e_u)-1] )
      if e_H[len(e_H)-1] > 2e-3:
         print('Python Warning: Pressure in shunn3 CC is out of tolerance. e_H = ',e_H[len(e_H)-1] )
