# General Plotting function.
# Each call to this module creates a single .pdf plot file.
# Input parameters required:
#   x_data - pre-processed data.
#   y_data - pre-processed data with mask_value set for missing data.
#   dictionary record for plot settings.

import matplotlib
matplotlib.use('PDF')
#matplotlib.use('SVG')
from matplotlib import rcParams
rcParams['font.family'] = 'serif'
rcParams['font.serif'] = ['Times New Roman']
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy.ma as M
from pylab import *
from math import sqrt

plotTypeVal = 's' # s = scatter, l = line
colorname = ['r','g','b']
symboltype = ['o','^','s']
linethickness = [1.0,0.5,1.0]
linetype = [':','-','-']
point_size = [30,20,80]
mask_value = -9999
xmin = 0
xmax = 100
ymin = 0
ymax = 100
output_filename = 'tester'
legend_loc='lower right'

if plotTypeVal == 's':
    x_data_set = [[40.00,60.00,125.00,165.00,245.00],[300.00, 430.00, 500.00, 555.00, 600.00],[620.00, 615.00, 600.00, 653.00, 660.00, 646.00]]
    y_data_set = [[50.80,69.23,119.71,167.59,256.20],[356.02, 425.82, 504.18, 553.45, 602.47],[624.12, 624.60, 656.12, 670.82, 668.81, 660.29]]
    xmin = 0
    xmax = 700
    ymin = 0
    ymax = 700
else:
    x_data_set = [[1.00,5.00,10.00,15.00,20.00,50.00,80.00],[10.00, 40.00, 50.00, 55.00, 60.00, 100.00],[1.00, 15.00, 60.00, 65.00, 70.00, 100.00]]
    y_data_set = [[50.80,69.23,-9999,167.59,256.20,450.00,600.00],[356.02, 425.82, -9999, 553.45, 602.47, 650.00],[624.12, 624.60, 656.12, 670.82, 668.81, 660.29]]
    xmin = 0
    xmax = 100
    ymin = 0
    ymax = 700

# Mask out missing or bad data in y_data_set
y_data_to_plot = []
for y_data_subset in y_data_set:
    y_values = M.array(y_data_subset)
    y_values_masked = M.masked_where(y_values == mask_value , y_values)
    y_data_to_plot.append(y_values_masked)

# Compute Scatter Plot Data Statistics
if plotTypeVal == 's':
    epsilion_vals = []
    for i in range(len(x_data_set)):
        for j in range(len(x_data_set[i])):
            if y_data_set[i][j] == mask_value:
                pass
            else:
                epsilion_vals = epsilion_vals+[(y_data_set[i][j]-x_data_set[i][j])/x_data_set[i][j]]
                
    mu_val = sum(epsilion_vals)/len(epsilion_vals)
    print 'Mu:', mu_val
    
    sigma_val = sqrt(sum([(e_val-mu_val)**2 for e_val in epsilion_vals])/len(epsilion_vals))
    sigma_2_val = sigma_val*2
    print '2 Sigma:', sigma_2_val
    
    mu_max = ymax*(1+mu_val)
    print 'Mu Max:', mu_max
#else:
#    pass
    
# Make Plot
if plotTypeVal == 's':
    fig = plt.figure(figsize=(5,5)) # (w,h) values are size in inches.
    ax = fig.add_subplot(111, aspect='equal')
else:
    fig = plt.figure(figsize=(6,4)) # (w,h) values are size in inches.
    ax = fig.add_subplot(111)

for pl in range(len(x_data_set)):
    if plotTypeVal == 's':
        ax.scatter(x_data_set[pl], y_data_to_plot[pl], s=point_size[pl], c=colorname[pl], marker=symboltype[pl], edgecolors='none', label='Data Set '+str(pl))
    else:
        ax.plot(x_data_set[pl], y_data_to_plot[pl], c=colorname[pl], linestyle=linetype[pl], linewidth=linethickness[pl], label='Data Set '+str(pl))

if plotTypeVal == 's':
    ax.plot([xmin,xmax],[ymin,ymax], 'k-', linewidth=0.5, label='_nolegend_') # Center Line
    muline = ax.plot([xmin,xmax],[ymin,mu_max], 'g:', linewidth=2.0, label='_nolegend_') # Mu Line
    sigma2upper = ax.plot([xmin,xmax], [ymin,(mu_max*(1+sigma_2_val))], 'r-.', linewidth=2.0, label='_nolegend_') # Upper 2 Sigma Line
    sigma2lower = ax.plot([xmin,xmax], [ymin,(mu_max*(1-sigma_2_val))], 'r-.', linewidth=2.0, label='_nolegend_') # Lower 2 Sigma Line
    plt.xlabel('Measured')
    plt.ylabel('Predicted')
    ax.text(xmax*0.015, ymax*0.94, 'Really Long Scatter Plot Title', horizontalalignment='left') 
    ax.text(xmax*0.030, ymax*0.90, r'$2\sigma='+"%2.2f" % (sigma_2_val*100)+'\%$')
    ax.text(xmax*0.030, ymax*0.85, r'$\mu='+"%2.2f" % (mu_val*100)+'\%$')
    # matplotlib.text.Text instances
    if legend_loc != '':
        leg = ax.legend(loc=legend_loc)
else:
    plt.xlabel('Time (s)')
    plt.ylabel('Data Values '+r'$^\circ$C')
    ax.text(xmax*0.015, ymax*0.94, 'Really Long Scatter Plot Title', horizontalalignment='left') 
    ax.grid(False)
    if legend_loc != '':
        leg = ax.legend(loc=legend_loc)

ax.axis([xmin, xmax, ymin, ymax])
if legend_loc != '':
    for t in leg.get_texts():
        t.set_fontsize('small')    # the legend text fontsize

plt.show()
plt.savefig(output_filename)

#return(0)