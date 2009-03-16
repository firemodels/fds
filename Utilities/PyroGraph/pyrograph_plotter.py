"""
pyrograph_plotter.py 

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
#rcParams['font.serif'] = ['Times New Roman']
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy.ma as M
from pylab import *
import pyrograph_calcs as calc
import pyrograph_parser as prsr


def scatter_plot(quantity_id,data_set,output_directory,diagnostic_level):
    # Read in quantities dictionary object.
    quantities = prsr.read_pickle("pyrograph_quantities_object.pkl",diagnostic_level)
    styles = prsr.read_pickle("pyrograph_styles_object.pkl",diagnostic_level)
    groups = prsr.read_pickle("pyrograph_groups_object.pkl",diagnostic_level)
    
    #Set Scatter Plot Settings
    title = quantities[quantity_id]["Scatter_Plot_Title"]
    size = float(quantities[quantity_id]["Plot_Width"])
    plot_min = float(quantities[quantity_id]["Plot_Min"])
    plot_max = float(quantities[quantity_id]["Plot_Max"])
    output_filename = quantities[quantity_id]["Plot_Filename"]
    legend_loc= quantities[quantity_id]["Key_Position"]
    ind_axis_title = quantities[quantity_id]["Ind_Title"]
    dep_axis_title = quantities[quantity_id]["Dep_Title"]
    title_position = eval(quantities[quantity_id]["Title_Position"])
    percent_error = float(quantities[quantity_id]["%error"])/100.0
    
    exp_error_lines = quantities[quantity_id]["Exp_Error_Lines"]
    mu_line = quantities[quantity_id]["Mu_Line"]
    sigma_lines = quantities[quantity_id]["Sigma_Lines"]
    
    #print "Plotting:",title
    
    x_data_set = []
    y_data_set = []
    
    fig = plt.figure(figsize=(float(size),float(size))) # (w,h) values are size in inches.
    ax = fig.add_subplot(111, aspect='equal')
    
    # Iterate through data_set and plot by Group_ID.
    for group_id in data_set:
        group_title = groups[group_id]["Group_Title"]
        symboltype = styles[int(groups[group_id]["Style_ID"])]['Symbol_Style']
        symbolcolor = styles[int(groups[group_id]["Style_ID"])]['Fill_Color']
        symbolsize = int(styles[int(groups[group_id]["Style_ID"])]['Symbol_Size'])
        edgecolor = styles[int(groups[group_id]["Style_ID"])]['Edge_Color']
        
        #print "Data:",data_set[group_id]
        x_data = [x[0] for x in data_set[group_id]]
        y_data = [x[1] for x in data_set[group_id]]
        x_data_set.append(x_data)
        y_data_set.append(y_data)
        ax.scatter(x_data, y_data, s=symbolsize, c=symbolcolor, marker=symboltype, edgecolors=edgecolor, label=group_title)
    
    mu_sigma = calc.mu_2sigma(x_data_set,y_data_set,diagnostic_level)
    mu_val = mu_sigma[0]
    sigma_2_val = mu_sigma[1]
    mu_max = plot_max*(1+mu_val)
    
    # Draw Center Line
    ax.plot([plot_min,plot_max],[plot_min,plot_max], 'k-', linewidth=2.0, label='_nolegend_')
    
    if exp_error_lines == 'yes':
        # Draw Upper Exp Error
        exp_error_upper = ax.plot([plot_min,plot_max], [plot_min,(plot_max*(1+percent_error))], 'k:', linewidth=1.0, label='_nolegend_')
        # Draw Lower Exp Error
        exp_error_lower = ax.plot([plot_min,plot_max], [plot_min,(plot_max*(1-percent_error))], 'k:', linewidth=1.0, label='_nolegend_')
    
    if mu_line == 'yes':
        # Draw Mu Line
        muline = ax.plot([plot_min,plot_max],[plot_min,mu_max], 'r-', linewidth=1.0, label='_nolegend_')
    
    if sigma_lines == 'yes':
        # Upper 2 Sigma Line
        sigma2upper = ax.plot([plot_min,plot_max], [plot_min,(mu_max*(1+sigma_2_val))], 'r:', linewidth=1.0, label='_nolegend_')
        # Lower 2 Sigma Line
        sigma2lower = ax.plot([plot_min,plot_max], [plot_min,(mu_max*(1-sigma_2_val))], 'r:', linewidth=1.0, label='_nolegend_')
    
    plt.xlabel(ind_axis_title)
    plt.ylabel(dep_axis_title)
    ax.text((plot_max-plot_min)*title_position[0], (plot_max-plot_min)*title_position[1], title, horizontalalignment='left') 
    ax.text((plot_max-plot_min)*(title_position[0]+0.05), (plot_max-plot_min)*(title_position[1]-0.07), r'$\mu='+"%2.2f" % (mu_val*100)+'\%$')
    ax.text((plot_max-plot_min)*(title_position[0]+0.05), (plot_max-plot_min)*(title_position[1]-0.12), r'$2\sigma='+"%2.2f" % (sigma_2_val*100)+'\%$')
    
    if legend_loc != '':
        leg = ax.legend(loc=legend_loc)
    ax.axis([plot_min, plot_max, plot_min, plot_max])
    if legend_loc != '':
        for t in leg.get_texts():
            t.set_fontsize('small')    # the legend text fontsize
            
    plt.show()
    plt.savefig(output_directory+output_filename)


def comparison_plot(data_set,data_info,d1_index_set,d2_index_set,output_directory,diagnostic_level):
    styles = prsr.read_pickle("pyrograph_styles_object.pkl",diagnostic_level)
    
    title = data_info["Plot_Title"]
    size = float(data_info["Plot_Width"])
    output_filename = data_info["Plot_Filename"]
    legend_loc= data_info["Key_Position"]
    ind_axis_title = data_info["Ind_Title"]
    dep_axis_title = data_info["Dep_Title"]
    
    xmin = float(data_info['Min_Ind'])
    xmax = float(data_info['Max_Ind'])
    ymin = float(data_info['Min_Dep'])
    ymax = float(data_info['Max_Dep'])
    
    title_position = eval(data_info["Title_Position"])
    title_x_pos = float(title_position[0])
    title_y_pos = float(title_position[1])
    
    mask_value = -9999.0
    
    if data_info['d1_Key'][0] == '[':
        d1_column_names = eval(data_info['d1_Key'])
    else:
        d1_column_names = [data_info['d1_Key']]
        
    if data_info['d2_Key'][0] == '[':
        d2_column_names = eval(data_info['d2_Key'])
    else:
        d2_column_names = [data_info['d2_Key']]
    
    fig = plt.figure(figsize=(size,size*0.75)) # (w,h) values are size in inches.
    ax = fig.add_subplot(111)
    
    plot_type_function = {'linear':ax.plot,'loglog':ax.loglog,'semilogx':ax.semilogx,'semilogy':ax.semilogy}
    
    d1_data = data_set[0]
    d2_data = data_set[1]
    d1_x_data = d1_data[0]
    d2_x_data = d2_data[0]
    
    d1_y_data = d1_data[1:]
    d2_y_data = d2_data[1:]
    
    # Mask out data set to -9999.0 in y_data_set
    d1_y_data_to_plot = []
    d2_y_data_to_plot = []
    
    for d1_y_data_subset in d1_y_data:
        y_values = M.array(d1_y_data_subset)
        d1_y_values_masked = M.masked_where(y_values == mask_value , y_values)
        d1_y_data_to_plot.append(d1_y_values_masked)
    
    for d2_y_data_subset in d2_y_data:
        y_values = M.array(d2_y_data_subset)
        d2_y_values_masked = M.masked_where(y_values == mask_value , y_values)
        d2_y_data_to_plot.append(d2_y_values_masked)
    
    linecount = 0
    for d1_y_masked_data in d1_y_data_to_plot:
        
        d1_style = eval(data_info["d1_Style"])
        d1_key = d1_column_names[linecount]
        symboltype = styles[int(d1_style[linecount])]['Symbol_Style']
        edgecolor = styles[int(d1_style[linecount])]['Edge_Color']
        symbolcolor = styles[int(d1_style[linecount])]['Fill_Color']
        
        if styles[int(d1_style[linecount])]['Symbol_Size'] == 'none':
            symbolsize = 'none'
        else:
            symbolsize = int(styles[int(d1_style[linecount])]['Symbol_Size'])
        
        linecolor = styles[int(d1_style[linecount])]['Line_Color']
        linestyle = styles[int(d1_style[linecount])]['Line_Style']
        
        if styles[int(d1_style[linecount])]['Line_Width'] == 'none':
            linewidth = 'none'
        else:
            linewidth = float(styles[int(d1_style[linecount])]['Line_Width'])
        
        plot_type_function[data_info['Plot_Type']](d1_x_data, d1_y_masked_data, c=linecolor, linestyle=linestyle, linewidth=linewidth, label=d1_key)
        linecount += 1
        
    linecount = 0
    for d2_y_masked_data in d2_y_data_to_plot:
        
        d2_style = eval(data_info["d2_Style"])
        d2_key = d2_column_names[linecount]
        symboltype = styles[int(d2_style[linecount])]['Symbol_Style']
        edgecolor = styles[int(d2_style[linecount])]['Edge_Color']
        symbolcolor = styles[int(d2_style[linecount])]['Fill_Color']
        
        if styles[int(d2_style[linecount])]['Symbol_Size'] == 'none':
            symbolsize = 'none'
        else:
            symbolsize = int(styles[int(d2_style[linecount])]['Symbol_Size'])
        
        linecolor = styles[int(d2_style[linecount])]['Line_Color']
        linestyle = styles[int(d2_style[linecount])]['Line_Style']
        
        if styles[int(d2_style[linecount])]['Line_Width'] == 'none':
            linewidth = 'none'
        else:
            linewidth = float(styles[int(d2_style[linecount])]['Line_Width'])
        
        plot_type_function[data_info['Plot_Type']](d2_x_data, d2_y_masked_data, c=linecolor, linestyle=linestyle, linewidth=linewidth, label=d2_key)
        linecount += 1
        
    plt.xlabel(ind_axis_title)
    plt.ylabel(dep_axis_title)
    
    ax.grid(False)
    ax.axis([xmin,xmax,ymin,ymax])
    
    ax.text((xmax-xmin)*title_x_pos, (ymax-ymin)*title_y_pos, title, horizontalalignment='left')
    
    if legend_loc != 'none':
        leg = ax.legend(loc=legend_loc)
    
    if legend_loc != 'none':
            for t in leg.get_texts():
                t.set_fontsize('small')    # the legend text fontsize
            
    plt.show()
    plt.savefig(output_directory+output_filename)

