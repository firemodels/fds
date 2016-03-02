#!/usr/bin/python
#McDermott
#2016-03-02

from __future__ import division # make floating point division default as in Matlab, e.g., 1/2=0.5
import math
import numpy as np
import scipy.special as sp
import matplotlib as mat
import matplotlib.pyplot as plt
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica'],'size':16})
mat.rcParams['mathtext.fontset'] = 'custom'
mat.rcParams['mathtext.rm'] = 'Helvetica'        #'Bitstream Vera Sans'
mat.rcParams['mathtext.it'] = 'Helvetica:italic' #'Bitstream Vera Sans:italic'
mat.rcParams['mathtext.bf'] = 'Helvetica:bold'   #'Bitstream Vera Sans:bold'


ddir = '../../Verification/Flowfields/'

M_16 = np.genfromtxt(ddir+'blasius_16_line.csv', delimiter=',', skip_header=1, names=True)
z_16 = M_16['Up1z']
u_16 = M_16['Up1']

M_32 = np.genfromtxt(ddir+'blasius_32_line.csv', delimiter=',', skip_header=1, names=True)
z_32 = M_32['Up1z']
u_32 = M_32['Up1']

M_64 = np.genfromtxt(ddir+'blasius_64_line.csv', delimiter=',', skip_header=1, names=True)
z_64 = M_64['Up1z']
u_64 = M_64['Up1']

# plot FDS results

plt.figure

marker_style = dict(color='red', linestyle=':', marker='o', fillstyle='none', markersize=0)
plt.plot(u_16,z_16,label='FDS $N$=16',**marker_style)

marker_style = dict(color='green', linestyle='-.', marker='o', fillstyle='none', markersize=0)
plt.plot(u_32,z_32,label='FDS $N$=32',**marker_style)

marker_style = dict(color='blue', linestyle='-', marker='o', fillstyle='none', markersize=0)
plt.plot(u_64,z_64,label='FDS $N$=64',**marker_style)

#plt.axis([min(t), max(t), min(HRR), 1.1*max(HRR)])
plt.xlabel('u (m/s)')
plt.ylabel('z (m)')
plt.legend(loc='upper left', numpoints=1, frameon=False)
#plt.show()
plt.savefig('blasius_prof.pdf', format='pdf')
plt.close()



# #gather blasius profile

# u0= max(u_64(:));
# zmax=0.3;
# [eta,fp] = blasius_analytic(u0, zmax);

# mu = 0.001;
# rho = 1.19987607036;
# xc = 0.05;
# z_blasius=eta*sqrt(mu/rho*xc/u0);
# u_blasius=fp*u0;


# %plot whole velocity profile
# range = 1:4:length(u_blasius);
# H(1)=plot(u_blasius(range),z_blasius(range),'bo');
# hold on
# H(2)=plot(u_16,z_16,'r--');
# H(3)=plot(u_32,z_32,'c-.');
# H(4)=plot(u_64,z_64,'g-');

# axis([0 1.1 0 0.15])

# set(gca,'Units',Plot_Units)
# set(gca,'FontName',Font_Name)
# set(gca,'FontSize',Title_Font_Size)
# set(gca,'Position',[Plot_X,Plot_Y,Plot_Width,Plot_Height])

# xlabel('{\it u} (m/s)','Interpreter',Font_Interpreter,'FontName',Font_Name,'FontSize',Label_Font_Size)
# ylabel('{\it z} (m)','Interpreter',Font_Interpreter,'FontName',Font_Name,'FontSize',Label_Font_Size)

# h = legend(H,'Blasius','{\it N_z}=16','{\it N_z}=32','{\it N_z}=64','Location','northwest');
# set(h,'FontSize',Title_Font_Size,'Interpreter',Font_Interpreter,'FontSize',Label_Font_Size)

# % add Git revision if file is available

# Git_Filename = [repository,'blasius_16_git.txt'];
# addverstr(gca,Git_Filename,'linear')

# % print to pdf for whole velocity profile
# set(gcf,'Visible',Figure_Visibility);
# set(gcf,'PaperUnits',Paper_Units);
# set(gcf,'PaperSize',[Paper_Width Paper_Height]);
# set(gcf,'PaperPosition',[0 0 Paper_Width Paper_Height]);
# print(gcf,'-dpdf','../../Manuals/FDS_Verification_Guide/SCRIPT_FIGURES/blasius_profile')


# %%%%%%%get error comparing with analytic solution(blasius)%%%%%%%%%%
# figure
# clear H

# err(1) = 0;
# err(2) = 0;
# err(3) = 0;

# % get error(n=16)
# for i=1:length(u_16)   
#     err(1)=err(1)+(abs(u_16(i)-u_blasius(9+(i-1)*16)))^2;    
# end    
# err(1)=err(1)/16;
# err(1)=sqrt(err(1));

# % get error(n=32)
# for i=1:length(u_32)    
#     err(2)=err(2)+(abs(u_32(i)-u_blasius(5+(i-1)*8)))^2 ;   
# end   
# err(2)=err(2)/32;
# err(2)=sqrt(err(2));

# % get error(n=64)
# for i=1:length(u_64)    
#     err(3)=err(3)+(abs(u_64(i)-u_blasius(3+(i-1)*4)))^2;    
# end   
# err(3)=err(3)/64;
# err(3)=sqrt(err(3));


# dz(1)=abs(z_16(10)-z_16(9));
# dz(2)=abs(z_32(10)-z_32(9));
# dz(3)=abs(z_64(10)-z_64(9));

# plot_style
# set(gca,'Units',Plot_Units)
# set(gca,'Position',[Plot_X,Plot_Y,Plot_Width,Plot_Height])

# H(1)=loglog(dz, err,'b*-','LineWidth',Line_Width); hold on
# H(2)=loglog(dz, 10*dz,'k--','LineWidth',Line_Width);
# H(3)=loglog(dz, 100*dz.^2,'k-','LineWidth',Line_Width);

# set(gca,'FontName',Font_Name)
# set(gca,'FontSize',Label_Font_Size)

# xlabel('Grid Spacing, {\it\deltaz} (m)','FontSize',Title_Font_Size,'Interpreter',Font_Interpreter,'FontSize',Label_Font_Size)
# ylabel('RMS Error (m/s)','FontSize',Title_Font_Size,'Interpreter',Font_Interpreter,'FontSize',Label_Font_Size)
# axis([1e-3 1e-1 1e-3 1e-0])
# legend_handle=legend(H,'FDS','{\itO}({\it\deltaz})','{\itO}({\it\deltaz^2})','Location','Northwest');
# set(legend_handle,'FontSize',Label_Font_Size,'Interpreter',Font_Interpreter)

# % add Git revision if file is available

# Git_Filename = [repository,'blasius_16_git.txt'];
# addverstr(gca,Git_Filename,'loglog')

# % print to pdf
# set(gcf,'Visible',Figure_Visibility);
# set(gcf,'PaperUnits',Paper_Units);
# set(gcf,'PaperSize',[Paper_Width Paper_Height]);
# set(gcf,'PaperPosition',[0 0 Paper_Width Paper_Height]);
# print(gcf,'-dpdf','../../Manuals/FDS_Verification_Guide/SCRIPT_FIGURES/blasius_convergence')



