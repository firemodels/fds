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
GIT=

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
rcParams.update({'font.size':20})
font1 = {'weight' : 'normal',
        'size'   : 22,
        }
font2 = {'weight' : 'normal',
        'size'   : 16,
        }

font3 = {'weight' : 'normal',
        'size'   : 14,
        }

# Defining adress for plots
fds_address='../examples/simple_panel_hot/'
fds_filename='simple_panel_hot'+'_devc.csv'
ansys_address='../examples/simple_panel_hot/'    
ansys_filename1='simple_panel_hot'+'.csv'
ansys_filename2='simple_panel_hot'+'_ast.csv'

# SIMPLE_PANEL_HOT CASE

# Plot 1 - simple_panel_hot - surface temperature

fds_result = pd.read_csv(fds_address+fds_filename, header=1)
ansys_result = pd.read_csv(ansys_address+ansys_filename1, header=0)

fig = figure()
plt.plot(fds_result['Time']/60,fds_result['T1'],'k-',mfc='none',label='FDS(T1)',linewidth=1)
plt.plot(fds_result['Time']/60,fds_result['T2'],'r-',mfc='none',label='FDS(T2)',linewidth=1)  
plt.plot(fds_result['Time']/60,fds_result['T3'],'b-',mfc='none',label='FDS(T3)',linewidth=1)    

plt.plot(ansys_result['Time']/60,ansys_result['T1'],'ko',mfc='none',mec='black',label='ANSYS(T1)',mew=1.15)
plt.plot(ansys_result['Time']/60,ansys_result['T2'],'ro',mfc='none',mec='red',label='ANSYS(T2)',mew=1.15)  
plt.plot(ansys_result['Time']/60,ansys_result['T3'],'bo',mfc='none',mec='blue',label='ANSYS(T3)',mew=1.15)    

plt.xlim([0, 10])
plt.ylim([0, 500])
plt.xlabel('Time (min)', fontdict=font1)
plt.ylabel('Temperature ($^{\circ}$C)', fontdict=font1)
plt.legend(numpoints=1,loc=4, prop=font2)
plt.text( 0.2, 465, 'Surface Temperature (simple_panel_hot)')
plt.text( 7.5, 507, GIT, fontdict=font3)
plt.savefig('SCRIPT_FIGURES/simple_panel_hot'+'_ts'+'.pdf',format='pdf')
close()

# calculating error
a=np.abs(ansys_result['T2'])
b=np.abs(fds_result['T2'])
result=1-a[19]/b[600]
tolerance=0.01
if result<=tolerance:
	within_tolerance='Yes'
else:
	within_tolerance='No'
results=([['simple\_panel\_hot',"{:1.3e}".format(b[600]),"{:1.3e}".format(a[19]),'Relative',"{:1.3e}".format(result),tolerance,within_tolerance]])


# Plot 2 - simple_panel_hot - adiabatic surface temperature
ansys_result = pd.read_csv(ansys_address+ansys_filename2, header=0)

fig = figure()
plt.plot(ansys_result['Time']/60,ansys_result['AST2'],'ro',mfc='none',mec='red',label='FTMI(T2)',mew=1.15)  

plt.plot(fds_result['Time']/60,fds_result['AST2'],'k-',mfc='none',label='FDS(T2)',linewidth=1)  

plt.xlim([0, 10])
plt.ylim([0, 750])
plt.xlabel('Time (min)', fontdict=font1)
plt.ylabel('Temperature ($^{\circ}$C)', fontdict=font1)
plt.legend(numpoints=1,loc=4, prop=font2)
plt.text( 0.4, 590, 'Adiabatic Surface Temperature (simple_panel_hot)')
plt.text( 7.5, 764, GIT, fontdict=font3)
plt.savefig('SCRIPT_FIGURES/simple_panel_hot'+'_ast'+'.pdf',format='pdf')
close()

# calculating error
x = np.abs(fds_result['Time'])
y = np.abs(fds_result['AST2'])
fds_ast=integrate.trapz(y,x)
x2 = np.abs(ansys_result['Time'])
y2 = np.abs(ansys_result['AST2'])
ansys_ast=integrate.trapz(y2,x2)
result=1-(ansys_ast/fds_ast)
tolerance=0.01
if result<=tolerance:
	within_tolerance='Yes'
else:
	within_tolerance='No'
results=np.concatenate((results,[['simple\_panel\_hot',"{:1.3e}".format(fds_ast),"{:1.3e}".format(ansys_ast),'Relative',"{:1.3e}".format(result),tolerance,within_tolerance]]))

# H_PROFILE CASE

# Plot 3 - h_profile - surface temperature
ansys_address='../examples/h_profile/'
ansys_filename1='h_profile'+'.csv'
ansys_result = pd.read_csv(ansys_address+ansys_filename1, header=0)

fig = figure()
plt.plot(ansys_result['Time']/60,ansys_result['hightemp_A'],'r-',mfc='none',label='Ts_A',linewidth=1)  
plt.plot(ansys_result['Time']/60,ansys_result['lowtemp_A'],'k--',mfc='none',label='Ts_A_back',linewidth=1)  
plt.plot(ansys_result['Time']/60,ansys_result['hightemp_B'],'b-',mfc='none',label='Ts_B',linewidth=1)  
plt.plot(ansys_result['Time']/60,ansys_result['lowtemp_B'],'k--',mfc='none',label='Ts_B_back',linewidth=1)  

plt.xlim([0, 60])
plt.ylim([0, 750])
plt.xlabel('Time (min)', fontdict=font1)
plt.ylabel('Temperature ($^{\circ}$C)', fontdict=font1)
plt.legend(numpoints=1,loc=4, prop=font2)
plt.text( 3, 40, 'Surface Temperature (h_profile)')
plt.text( 45, 760, GIT, fontdict=font3)
plt.savefig('SCRIPT_FIGURES/h_profile'+'_ts'+'.pdf',format='pdf')
close()

# calculating error
a=np.abs(ansys_result['hightemp_B'])[120]
b=np.abs(ansys_result['lowtemp_B'])[120]
result=a-b
tolerance=0
if result>=tolerance:
	within_tolerance='Yes'
else:
	within_tolerance='No'
results=np.concatenate((results,[['h\_profile','> '+"{:1.3e}".format(b),"{:1.3e}".format(a),'Logic',' - ',' - ',within_tolerance]]))


# Plot 4 - h_profile - horizontal displacement (x)
fig = figure()
plt.plot(ansys_result['Time']/60,ansys_result['dx_A']*100,'r-',mfc='none',label='dx_A',linewidth=1.5)  
plt.plot(ansys_result['Time']/60,ansys_result['dx_B']*100,'k-',mfc='none',label='dx_B',linewidth=1.5)  
plt.plot(ansys_result['Time']/60,ansys_result['dx_C']*100,'b--',mfc='none',label='dx_C',linewidth=1.5)  

plt.xlim([0, 60])
plt.ylim([0, 2])
plt.xlabel('Time (min)', fontdict=font1)
plt.ylabel('Displacement (cm)', fontdict=font1)
plt.legend(numpoints=1,loc=4, prop=font2)
plt.text( 2, 1.83, 'Displacement in x axis (h_profile)')
plt.text( 45, 2.02, GIT, fontdict=font3)
plt.savefig('SCRIPT_FIGURES/h_profile'+'_dx'+'.pdf',format='pdf')
close()

# calculating error
a = np.abs(ansys_result['dx_A'])
b = np.abs(ansys_result['dx_B'])
result=np.average(a-b)
tolerance=0.01
if result<=tolerance:
	within_tolerance='Yes'
else:
	within_tolerance='No'
results=np.concatenate((results,[['h\_profile',"{:1.3e}".format(np.average(a)),"{:1.3e}".format(np.average(b)),'Absolute',"{:1.3e}".format(result),tolerance,within_tolerance]]))

# Plot 5 - h_profile - horizontal displacement (y)
fig = figure()
plt.plot(ansys_result['Time']/60,ansys_result['dy_A']*100,'r-',mfc='none',label='dy_A',linewidth=1.5)  
plt.plot(ansys_result['Time']/60,ansys_result['dy_C']*100,'b-',mfc='none',label='dy_C',linewidth=1.5)  

plt.xlim([0, 60])
plt.ylim([-0.3, 0.3])
plt.xlabel('Time (min)', fontdict=font1)
plt.ylabel('Displacement (cm)', fontdict=font1)
plt.legend(numpoints=1,loc=7, prop=font2)
plt.text( 5, -0.01, 'Displacement in y axis (h_profile)')
plt.text( 45, 0.308, GIT, fontdict=font3)
plt.savefig('SCRIPT_FIGURES/h_profile'+'_dy'+'.pdf',format='pdf')
close()

# calculating error
a = np.abs(ansys_result['dy_A'])
b = np.abs(ansys_result['dy_C'])
result=np.average(a+b)
tolerance=0.01
if result<=tolerance:
	within_tolerance='Yes'
else:
	within_tolerance='No'
results=np.concatenate((results,[['h\_profile',"{:1.3e}".format(np.average(a)),"{:1.3e}".format(np.average(b)),'Absolute',"{:1.3e}".format(result),tolerance,within_tolerance]]))
np.savetxt("verification_statistics.tex", results, delimiter="&", fmt="%s", newline='\\\\ \n')