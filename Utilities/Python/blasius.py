#!/usr/bin/python
#McDermott
#2016-03-02

from __future__ import division # make floating point division default as in Matlab, e.g., 1/2=0.5
import blaslib as bl
import math
import numpy as np
import scipy.special as sp
from scipy.interpolate import interp1d
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica'],'size':16})
matplotlib.rcParams['mathtext.fontset'] = 'custom'
matplotlib.rcParams['mathtext.rm'] = 'Helvetica'        #'Bitstream Vera Sans'
matplotlib.rcParams['mathtext.it'] = 'Helvetica:italic' #'Bitstream Vera Sans:italic'
matplotlib.rcParams['mathtext.bf'] = 'Helvetica:bold'   #'Bitstream Vera Sans:bold'


ddir = '../../Verification/Flowfields/'

M_16 = np.genfromtxt(ddir+'blasius_16_line.csv', delimiter=',', skip_header=1, names=True)
z_16 = M_16['Up05z']
u_16 = M_16['Up05']

M_32 = np.genfromtxt(ddir+'blasius_32_line.csv', delimiter=',', skip_header=1, names=True)
z_32 = M_32['Up05z']
u_32 = M_32['Up05']

# M_64 = np.genfromtxt(ddir+'blasius_64_line.csv', delimiter=',', skip_header=1, names=True)
# z_64 = M_64['Up05z']
# u_64 = M_64['Up05']

u0  = u_32[-1]
mu  = 0.001
rho = 1.199
xc  = 0.05

delta = math.sqrt(mu/rho*xc/u0)
zmax = 0.3
etamax = zmax/delta
n = 256
deta = etamax/n

[eta, fp] = bl.blasius_analytic(u0,zmax,mu,rho,xc)

z_bl=eta*delta
u_bl=fp*u0

# create a cubic interpolation from samples of the Blasius exact solution
# to be evaluated at FDS cell face centroid locations, where the u velocity
# component is stored.

u_interp = interp1d(z_bl, u_bl, kind='cubic')

# plot FDS results

plt.figure

marker_style = dict(color='red', linestyle='--', marker='o', fillstyle='none', markersize=5)
plt.plot(u_16,z_16,label='FDS $N$=16',**marker_style)

marker_style = dict(color='green', linestyle='--', marker='+', fillstyle='none', markersize=5)
plt.plot(u_32,z_32,label='FDS $N$=32',**marker_style)

# marker_style = dict(color='blue', linestyle='--', marker='^', fillstyle='none', markersize=5)
# plt.plot(u_64,z_64,label='FDS $N$=64',**marker_style)

# here is the exact Blasius solution
marker_style = dict(color='black', linestyle='-', marker='o', fillstyle='none', markersize=0)
plt.plot(u_bl,z_bl,label='Blasius exact',**marker_style)

# here is the Blasius solution sampled at the location of the N=16 u velocity component locations
marker_style = dict(color='red', linestyle='', marker='o', fillstyle='none', markersize=5)
plt.plot(u_interp(z_16),z_16,label='Blasius sample 16',**marker_style)

# here is the Blasius solution sampled at the location of the N=32 u velocity component locations
marker_style = dict(color='green', linestyle='', marker='+', fillstyle='none', markersize=5)
plt.plot(u_interp(z_32),z_32,label='Blasius sample 32',**marker_style)

plt.axis([0., 1.1, 0., 0.10])
plt.xlabel('u (m/s)')
plt.ylabel('z (m)')
plt.legend(loc='upper left', numpoints=1, frameon=False)
#plt.show()
plt.savefig('blasius_prof.pdf', format='pdf')
plt.close()

# error vectors

error_16 = u_16 - u_interp(z_16) # this creates a vector of errors of length len(z_16)
error_32 = u_32 - u_interp(z_32)

# RMS error

rms_16 = math.sqrt(np.linalg.norm(error_16)/len(z_16))
rms_32 = math.sqrt(np.linalg.norm(error_32)/len(z_32))

print rms_16, rms_32









