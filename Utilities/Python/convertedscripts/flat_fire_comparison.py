import numpy as np
import re
import smop_util
import matplotlib.pyplot as plt
# Trettel
# flat fire comparison
# This script just computes the exact solution. The actual plots are made by dataplot, and the exact solution is committed to the repo.

plt.close('all')

V_0=400
# temp/flat_fire_comparison.m:8
h=8
# temp/flat_fire_comparison.m:9
g=9.8
# temp/flat_fire_comparison.m:10
C_d=0.2
# temp/flat_fire_comparison.m:11
rho_a=1.2
# temp/flat_fire_comparison.m:12
rho_d=1000
# temp/flat_fire_comparison.m:13
D=0.005
# temp/flat_fire_comparison.m:14
K=np.dot(np.dot(3,rho_a),C_d) / (np.dot(np.dot(4,rho_d),D))
# temp/flat_fire_comparison.m:16
dt=0.05
# temp/flat_fire_comparison.m:18
tend=1.65
# temp/flat_fire_comparison.m:19
tvec=np.arange(0,tend+dt,dt)
# temp/flat_fire_comparison.m:21
xexact=np.log(np.dot(np.dot(V_0,K),tvec) + 1) / K
# temp/flat_fire_comparison.m:23
yexact=h + np.dot((g / (np.dot(2,(np.dot(K,V_0)) ** 2))),np.log(np.dot(np.dot(V_0,K),tvec) + 1)) - np.dot(g,tvec ** 2) / 4 - np.dot(g,tvec) / (np.dot(np.dot(2,K),V_0))
# temp/flat_fire_comparison.m:24
uexact=V_0 / (np.dot(np.dot(V_0,K),tvec) + 1)
# temp/flat_fire_comparison.m:28
vexact=g / (np.dot(np.dot(np.dot(2,K),V_0),(np.dot(np.dot(K,V_0),tvec) + 1))) - np.dot(g,tvec) / 2 - g / (np.dot(np.dot(2,K),V_0))
# temp/flat_fire_comparison.m:29
mat=np.asarray([], dtype='object')
# temp/flat_fire_comparison.m:33
mat = smop_util.safe_set(mat,('shape',0,),tvec)
# temp/flat_fire_comparison.m:34
mat = smop_util.safe_set(mat,('shape',1,),xexact)
# temp/flat_fire_comparison.m:35
mat = smop_util.safe_set(mat,('shape',2,),yexact)
# temp/flat_fire_comparison.m:36
mat = smop_util.safe_set(mat,('shape',3,),uexact)
# temp/flat_fire_comparison.m:37
mat = smop_util.safe_set(mat,('shape',4,),vexact)
# temp/flat_fire_comparison.m:38
np.savetxt('flat_fire.csv',mat,delimiter=',')