# -*- coding: utf-8 -*-
"""
Created on Sat Apr 23 13:18:18 2022

@author: Julio Cesar Silva
"""
from __future__ import division

import numpy as np
import pandas as pd
import matplotlib
matplotlib.get_configdir()
matplotlib.use('Agg') # Must be before importing matplotlib.pyplot or pylab!
import matplotlib.pyplot as plt

from pylab import *
from matplotlib import rcParams
from scipy import integrate

# Getting SVN number 
GIT='verification'

# Help to find the fonts properly!
# You can find the list of fonts installed on the system with:
#fonts = "\n".join(matplotlib.font_manager.findSystemFonts(fontpaths=None, fontext="ttf"))
#print fonts
# If Times New Roman is listed in the system, check if it is listed in python:
#print [f.name for f in matplotlib.font_manager.fontManager.ttflist]
# Maybe the python package have its own fonts folder, find it and add the file times.ttf
# Then, you need to remove the fontmanager cache so the code will read the new file
# Find the cache location with:
#print matplotlib.get_cachedir()

# Defining fonts to be used in plots
rcParams.update({'figure.autolayout': True})
rcParams.update({'font.family':'Times New Roman'})
rcParams.update({'font.size':16})
font1 = {'weight' : 'normal',
        'size'   : 18,
        }
font2 = {'weight' : 'normal',
        'size'   : 12,
        }

font3 = {'weight' : 'normal',
        'size'   : 10,
        }

# Defining adress for plots
fds_address=''
fds_filename='simple_panel'+'_devc.csv'
ftmi_address=''    

# SIMPLE_PANEL_HOT CASE

fds_result = pd.read_csv(fds_address+fds_filename, header=1)

# READING FTMI RESULT FROM DAT FILE
datfile=open(ftmi_address+'simple_panel_verification.dat', 'r')
ftmi_result=pd.DataFrame(columns=['Time','AST1','AST2','AST3'])  
i=0
j=0
line=0
while i <= 72:
    a=datfile.readline()
    i=i+1
while i <= 673:
    ftmi_result.at[j,'Time']=float(datfile.readline()[27:40])
    ftmi_result.at[j,'AST1']=float(datfile.readline()[27:40])
    j=j+1
    i=i+2
j=0
while i <= 1276:
    a=datfile.readline()
    i=i+1
while i <= 1877:
    ftmi_result.at[j,'Time']=float(datfile.readline()[27:40])
    ftmi_result.at[j,'AST2']=float(datfile.readline()[27:40])
    j=j+1   
    i=i+2
j=0
while i <= 2480:
    a=datfile.readline()
    i=i+1
while i <= 3081:
    ftmi_result.at[j,'Time']=float(datfile.readline()[27:40])
    ftmi_result.at[j,'AST3']=float(datfile.readline()[27:40])
    j=j+1   
    i=i+2    


# Plot simple_panel - adiabatic surface temperature

fig = plt.figure()
plt.plot(ftmi_result['Time']/60,ftmi_result['AST2'],'ro',mfc='none',mec='red',label='FTMI(AST2)',mew=1.15)  
plt.plot(fds_result['Time']/60,fds_result['AST2'],'k-',mfc='none',label='FDS(AST2)',linewidth=1)  
#plt.plot(ftmi_result['Time']/60,ftmi_result['AST3'],'ro',mfc='none',mec='red',label='FTMI(AST3)',mew=1.15)  
#plt.plot(fds_result['Time']/60,fds_result['AST3'],'k-',mfc='none',label='FDS(AST3)',linewidth=1)  
#plt.plot(ftmi_result['Time']/60,ftmi_result['AST1'],'ro',mfc='none',mec='red',label='FTMI(AST1)',mew=1.15)  
#plt.plot(fds_result['Time']/60,fds_result['AST1'],'k-',mfc='none',label='FDS(AST1)',linewidth=1)  

plt.xlim([0, 10])
plt.ylim([0, 750])
plt.xlabel('Time (min)', fontdict=font1)
plt.ylabel('Temperature ($^{\circ}$C)', fontdict=font1)
plt.legend(numpoints=1,loc=4, prop=font2)
plt.text( 0.4, 590, 'Adiabatic Surface Temperature (simple_panel)')
plt.text( 7.5, 764, GIT, fontdict=font3)
plt.savefig('simple_panel'+'_ast'+'.pdf',format='pdf')
plt.close()

# calculating error
x = np.abs(fds_result['Time'])
y1 = np.abs(fds_result['AST1'])
y2 = np.abs(fds_result['AST2'])
y3 = np.abs(fds_result['AST3'])
fds_ast1=integrate.trapz(y1,x)
fds_ast2=integrate.trapz(y2,x)
fds_ast3=integrate.trapz(y3,x)

x1 = np.abs(ftmi_result['Time'])
y11 = np.abs(ftmi_result['AST1'])
y12 = np.abs(ftmi_result['AST2'])
y13 = np.abs(ftmi_result['AST3'])
ftmi_ast1=integrate.trapz(y11,x1)
ftmi_ast2=integrate.trapz(y12,x1)
ftmi_ast3=integrate.trapz(y13,x1)

result=max(abs(1-(ftmi_ast1/fds_ast1)),abs(1-(ftmi_ast2/fds_ast2)),abs(1-(ftmi_ast3/fds_ast3)))
tolerance=0.025
if abs(result)<=tolerance:
	within_tolerance='Yes'
else:
	within_tolerance='No'
    
results=([['simple\_panel',"{:1.3e}".format(fds_ast2),"{:1.3e}".format(ftmi_ast2),'Relative',"{:1.3e}".format(result),tolerance,within_tolerance]])

# Plot simple_panel - convective heat transfer coefficient

# READING FTMI RESULT FROM DAT FILE
datfile=open(ftmi_address+'simple_panel_verification.dat', 'r')
ftmi_result=pd.DataFrame(columns=['Time','H1','H2','H3'])  
i=0
j=0
line=0
while i <= 674:
    a=datfile.readline()
    i=i+1
while i <= 1275:
    ftmi_result.at[j,'Time']=float(datfile.readline()[27:40])
    ftmi_result.at[j,'H1']=float(datfile.readline()[27:40])
    j=j+1
    i=i+2
j=0
while i <= 1878:
    a=datfile.readline()
    i=i+1
while i <= 2479:
    ftmi_result.at[j,'Time']=float(datfile.readline()[27:40])
    ftmi_result.at[j,'H2']=float(datfile.readline()[27:40])
    j=j+1   
    i=i+2
j=0
while i <= 3082:
    a=datfile.readline()
    i=i+1
while i <= 3683:
    ftmi_result.at[j,'Time']=float(datfile.readline()[27:40])
    ftmi_result.at[j,'H3']=float(datfile.readline()[27:40])
    j=j+1   
    i=i+2    


fig = plt.figure()
plt.plot(ftmi_result['Time']/60,ftmi_result['H2'],'ro',mfc='none',mec='red',label='FTMI(H2)',mew=1.15)  
plt.plot(fds_result['Time']/60,fds_result['H2'],'k-',mfc='none',label='FDS(H2)',linewidth=1)  
#plt.plot(ftmi_result['Time']/60,ftmi_result['H1'],'ro',mfc='none',mec='red',label='FTMI(H1)',mew=1.15)  
#plt.plot(fds_result['Time']/60,fds_result['H1'],'k-',mfc='none',label='FDS(H1)',linewidth=1)  
#plt.plot(ftmi_result['Time']/60,ftmi_result['H3'],'ro',mfc='none',mec='red',label='FTMI(H3)',mew=1.15)  
#plt.plot(fds_result['Time']/60,fds_result['H3'],'k-',mfc='none',label='FDS(H3)',linewidth=1)  

plt.xlim([0, 10])
plt.ylim(bottom=0)
plt.xlabel('Time (min)', fontdict=font1)
plt.ylabel('h (W/(m²K))', fontdict=font1)
plt.legend(numpoints=1,loc=4, prop=font2)
plt.text( 0.4, 590, 'Convective Heat Transfer (simple_panel)')
plt.text( 7.5, 764, GIT, fontdict=font3)
plt.savefig('simple_panel'+'_h'+'.pdf',format='pdf')
plt.close()

# calculating error
x = np.abs(fds_result['Time'])
y1 = np.abs(fds_result['H1'])
y2 = np.abs(fds_result['H2'])
y3 = np.abs(fds_result['H3'])
fds_h1=integrate.trapz(y1,x)
fds_h2=integrate.trapz(y2,x)
fds_h3=integrate.trapz(y3,x)

x1 = np.abs(ftmi_result['Time'])
y11 = np.abs(ftmi_result['H1'])
y12 = np.abs(ftmi_result['H2'])
y13 = np.abs(ftmi_result['H3'])
ftmi_h1=integrate.trapz(y11,x1)
ftmi_h2=integrate.trapz(y12,x1)
ftmi_h3=integrate.trapz(y13,x1)

result=max(abs(1-(ftmi_h1/fds_h1)),abs(1-(ftmi_h2/fds_h2)),abs(1-(ftmi_h3/fds_h3)))
tolerance=0.05
if abs(result)<=tolerance:
	within_tolerance='Yes'
else:
	within_tolerance='No'

results=np.concatenate((results,[['simple\_panel',"{:1.3e}".format(fds_h1),"{:1.3e}".format(ftmi_h1),'Relative',"{:1.3e}".format(result),tolerance,within_tolerance]]))
np.savetxt("verification_statistics.tex", results, delimiter="&", fmt="%s", newline='\\\\ \n')
