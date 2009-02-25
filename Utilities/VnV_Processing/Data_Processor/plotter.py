"""
General Plotting functions.

Input parameters required:
    x_data - pre-processed data.
    y_data - pre-processed data with mask_value set for missing data.
    dictionary record for plot settings.

"""

import matplotlib
matplotlib.use('PDF')
from matplotlib import rcParams
rcParams['font.family'] = 'serif'
rcParams['font.serif'] = ['Times New Roman']
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy.ma as M
from pylab import *
import calcs as clc
import parser as pr

def scatter_plot(quantity_id,data_set,output_directory,diagnostic_level):
    # Read in quantities dictionary object.
    quantities = pr.read_pickle("quantities_object.pkl",diagnostic_level)
    styles = pr.read_pickle("styles_object.pkl",diagnostic_level)
    groups = pr.read_pickle("groups_object.pkl",diagnostic_level)
    
    #Set Scatter Plot Settings
    title = quantities[quantity_id]["Scatter_Plot_Title"]
    size = float(quantities[quantity_id]["Plot_Width"])
    plot_min = float(quantities[quantity_id]["Plot_Min"])
    plot_max = float(quantities[quantity_id]["Plot_Max"])
    output_filename = quantities[quantity_id]["Plot_Filename"]
    legend_loc= quantities[quantity_id]["Key_Position"]
    ind_axis_title = quantities[quantity_id]["Ind_Title"]
    dep_axis_title = quantities[quantity_id]["Ind_Title"]
    title_position = eval(quantities[quantity_id]["Title_Position"])
    percent_error = int(quantities[quantity_id]["%error"])
    mask_value = -9999
    
    print "Plot Title:",title
    
    x_data_set = []
    y_data_set = []
    
    fig = plt.figure(figsize=(float(size),float(size))) # (w,h) values are size in inches.
    ax = fig.add_subplot(111, aspect='equal')
    
    # Iterate through data_set and plot by Group_ID.
    for group_id in data_set:
        #print "Group #:",group_id
        group_title = groups[group_id]["Group_Title"]
        symboltype = styles[int(groups[group_id]["Style_ID"])]['Symbol_Style']
        colorname = styles[int(groups[group_id]["Style_ID"])]['Fill_Color']
        symbolsize = int(styles[int(groups[group_id]["Style_ID"])]['Symbol_Size'])
        edgecolor = styles[int(groups[group_id]["Style_ID"])]['Edge_Color']
        
        #print "Data:",data_set[group_id]
        x_data = [x[0] for x in data_set[group_id]]
        y_data = [x[1] for x in data_set[group_id]]
        x_data_set.append(x_data)
        y_data_set.append(y_data)
        ax.scatter(x_data, y_data, s=symbolsize, c=colorname, marker=symboltype, edgecolors=edgecolor, label=group_title)

    
    mu_sigma = clc.mu_2sigma(x_data_set,y_data_set)
    mu_val = mu_sigma[0]
    sigma_2_val = mu_sigma[1]
    mu_max = plot_max*(1+mu_val)
    
    # Draw Center Line
    ax.plot([plot_min,plot_max],[plot_min,plot_max], 'k-', linewidth=0.5, label='_nolegend_')
    # Draw Mu Line
    muline = ax.plot([plot_min,plot_max],[plot_min,mu_max], 'g:', linewidth=2.0, label='_nolegend_')
    # Upper 2 Sigma Line
    sigma2upper = ax.plot([plot_min,plot_max], [plot_min,(mu_max*(1+sigma_2_val))], 'r-.', linewidth=2.0, label='_nolegend_')
    # Lower 2 Sigma Line
    sigma2lower = ax.plot([plot_min,plot_max], [plot_min,(mu_max*(1-sigma_2_val))], 'r-.', linewidth=2.0, label='_nolegend_')
    plt.xlabel(ind_axis_title)
    plt.ylabel(dep_axis_title)
    ax.text(plot_max*title_position[0], plot_max*title_position[1], title, horizontalalignment='left') 
    ax.text(plot_max*(title_position[0]+0.05), plot_max*(title_position[1]-0.05), r'$2\sigma='+"%2.2f" % (sigma_2_val*100)+'\%$')
    ax.text(plot_max*(title_position[0]+0.05), plot_max*(title_position[1]-0.10), r'$\mu='+"%2.2f" % (mu_val*100)+'\%$')
    
    if legend_loc != '':
        leg = ax.legend(loc=legend_loc)
    ax.axis([plot_min, plot_max, plot_min, plot_max])
    if legend_loc != '':
        for t in leg.get_texts():
            t.set_fontsize('small')    # the legend text fontsize
            
    plt.show()
    plt.savefig(output_directory+output_filename)
# 
# def comparison_plot():
#     
# 
# for pl in range(len(x_data_set)):
#     if plotTypeVal == 's':
#         ax.scatter(x_data_set[pl], y_data_to_plot[pl], s=point_size[pl], c=colorname[pl], marker=symboltype[pl], edgecolors='none', label='Data Set '+str(pl))
#     else:
#         ax.plot(x_data_set[pl], y_data_to_plot[pl], c=colorname[pl], linestyle=linetype[pl], linewidth=linethickness[pl], label='Data Set '+str(pl))
# 
# if plotTypeVal == 's':
#     mu_sigma = clc.mu2sigma()
#     mu_val = mu_sigma[0]
#     sigma_2_val = mu_sigma[1]
#     mu_max = 
#     
#     ax.plot([xmin,xmax],[ymin,ymax], 'k-', linewidth=0.5, label='_nolegend_') # Center Line
#     muline = ax.plot([xmin,xmax],[ymin,mu_max], 'g:', linewidth=2.0, label='_nolegend_') # Mu Line
#     sigma2upper = ax.plot([xmin,xmax], [ymin,(mu_max*(1+sigma_2_val))], 'r-.', linewidth=2.0, label='_nolegend_') # Upper 2 Sigma Line
#     sigma2lower = ax.plot([xmin,xmax], [ymin,(mu_max*(1-sigma_2_val))], 'r-.', linewidth=2.0, label='_nolegend_') # Lower 2 Sigma Line
#     plt.xlabel('Measured')
#     plt.ylabel('Predicted')
#     ax.text(xmax*0.015, ymax*0.94, 'Really Long Scatter Plot Title', horizontalalignment='left') 
#     ax.text(xmax*0.030, ymax*0.90, r'$2\sigma='+"%2.2f" % (sigma_2_val*100)+'\%$')
#     ax.text(xmax*0.030, ymax*0.85, r'$\mu='+"%2.2f" % (mu_val*100)+'\%$')
#     # matplotlib.text.Text instances
#     if legend_loc != '':
#         leg = ax.legend(loc=legend_loc)
# else:
#     plt.xlabel('Time (s)')
#     plt.ylabel('Data Values '+r'$^\circ$C')
#     ax.text(xmax*0.015, ymax*0.94, 'Really Long Scatter Plot Title', horizontalalignment='left') 
#     ax.grid(False)
#     if legend_loc != '':
#         leg = ax.legend(loc=legend_loc)
# 
# ax.axis([xmin, xmax, ymin, ymax])
# if legend_loc != '':
#     for t in leg.get_texts():
#         t.set_fontsize('small')    # the legend text fontsize
# 
# plt.show()
# plt.savefig(output_filename)
# 
# def test_scatter_plot():
#     colorname = ['r','g','b']
#     symboltype = ['o','^','s']
#     linethickness = [1.0,0.5,1.0]
#     linetype = [':','-','-']
#     point_size = [30,20,80]
#     mask_value = -9999
#     output_filename = 'tester'
#     legend_loc='lower right'
#     x_data_set = [[40.00,60.00,125.00,165.00,245.00],[300.00, 430.00, 500.00, 555.00, 600.00],[620.00, 615.00, 600.00, 653.00, 660.00, 646.00]]
#     y_data_set = [[50.80,69.23,119.71,167.59,256.20],[356.02, 425.82, 504.18, 553.45, 602.47],[624.12, 624.60, 656.12, 670.82, 668.81, 660.29]]
#     xmin = 0
#     xmax = 700
#     ymin = 0
#     ymax = 700
#     fig = plt.figure(figsize=(5,5)) # (w,h) values are size in inches.
#     ax = fig.add_subplot(111, aspect='equal')
#     
# def test_comparison_plot():
#     colorname = ['r','g','b']
#     symboltype = ['o','^','s']
#     linethickness = [1.0,0.5,1.0]
#     linetype = [':','-','-']
#     point_size = [30,20,80]
#     mask_value = -9999
#     output_filename = 'tester'
#     legend_loc='lower right'
#     x_data_set = [[1.00,5.00,10.00,15.00,20.00,50.00,80.00],[10.00, 40.00, 50.00, 55.00, 60.00, 100.00],[1.00, 15.00, 60.00, 65.00, 70.00, 100.00]]
#     y_data_set = [[50.80,69.23,-9999,167.59,256.20,450.00,600.00],[356.02, 425.82, -9999, 553.45, 602.47, 650.00],[624.12, 624.60, 656.12, 670.82, 668.81, 660.29]]
#     xmin = 0
#     xmax = 100
#     ymin = 0
#     ymax = 700
#     
#     fig = plt.figure(figsize=(6,4)) # (w,h) values are size in inches.
#     ax = fig.add_subplot(111)
#     
#     # Mask out missing or bad data in y_data_set
#     y_data_to_plot = []
#     for y_data_subset in y_data_set:
#         y_values = M.array(y_data_subset)
#         y_values_masked = M.masked_where(y_values == mask_value , y_values)
#         y_data_to_plot.append(y_values_masked)
#     
#     
# #return(0)