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
fds_address='../examples/simply_beam/'
#fds_filename='simply_beam'+'_devc.csv'
ansys_address='../examples/simply_beam/'    
ansys_filename1='simply_beam'+'.csv'


# SIMPLY_BEAM CASE

# Plot 1 - surface temperature at midspam and 0.5m from midspam

#fds_result = pd.read_csv(fds_address+fds_filename, header=1)
ansys_result = pd.read_csv(ansys_address+ansys_filename1, header=0)

fig = figure()
plt.plot(ansys_result['Time']/60,ansys_result['Temp_A'],'ko',mfc='none',mec='black',label='ANSYS(Temp_A)',mew=1.15)
plt.plot(ansys_result['Time']/60,ansys_result['Temp_B'],'ro',mfc='none',mec='red',label='ANSYS(Temp_B)',mew=1.15)  

plt.xlim([0, 10])
plt.ylim([0, 700])
plt.xlabel('Time (min)', fontdict=font1)
plt.ylabel('Temperature ($^{\circ}$C)', fontdict=font1)
plt.legend(numpoints=1,loc=4, prop=font2)
plt.text( 0.2, 650, 'Surface Temperature (simply_beam)')
plt.text( 6.8, 705, GIT, fontdict=font3)
plt.savefig('SCRIPT_FIGURES/simply_beam'+'_ts'+'.pdf',format='pdf')
close()


# Plot 2 - simply_beam - displacements in z

fig = figure()
plt.plot(ansys_result['Time']/60,ansys_result['dz_A']*100,'ko',mfc='none',mec='black',label='ANSYS(dz_A)',mew=1.15)  

plt.xlim([0, 10])
plt.ylim([-5, 0])
plt.xlabel('Time (min)', fontdict=font1)
plt.ylabel('Displacement (cm)', fontdict=font1)
plt.legend(numpoints=1,loc=4, prop=font2)
plt.text( 3.5, -0.35, 'Displacement in z axis (simply_beam)')
plt.text( 6.8, 0.05, GIT, fontdict=font3)
plt.savefig('SCRIPT_FIGURES/simply_beam'+'_dz'+'.pdf',format='pdf')
close()


# Plot 3 - simply_beam - displacements in x and y

fig = figure()
plt.plot(ansys_result['Time']/60,ansys_result['dx_A']*100,'ko',mfc='none',mec='black',label='ANSYS(dx_A)',mew=1.15)  
plt.plot(ansys_result['Time']/60,ansys_result['dy_A']*100,'ro',mfc='none',mec='red',label='ANSYS(dy_A)',mew=1.15)  
plt.xlim([0, 10])
plt.ylim([-0.1, 0.5])
plt.xlabel('Time (min)', fontdict=font1)
plt.ylabel('Displacement (cm)', fontdict=font1)
plt.legend(numpoints=1,loc=5, prop=font2)
plt.text( 0.2, 0.455, 'Displacement in x and y axis (simply_beam)')
plt.text( 6.8, 0.505, GIT, fontdict=font3)
plt.savefig('SCRIPT_FIGURES/simply_beam'+'_dx_dy'+'.pdf',format='pdf')
close()

