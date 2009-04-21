"""
pyrograph_plotter.py 

"""

import matplotlib
from matplotlib import rcParams
rcParams['font.family'] = 'serif'
rcParams['font.serif'] = ['Times New Roman']
rcParams['legend.fancybox'] = True
rcParams['legend.labelspacing'] = 0.2
#rcParams['pdf.use14corefonts'] = True
rcParams['font.size'] = 11
rcParams['axes.labelsize'] = 11
rcParams['legend.fontsize'] = 8
rcParams['xtick.labelsize'] = 11
rcParams['ytick.labelsize'] = 11

matplotlib.use('PDF')
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy.ma as M
from pylab import *
import modules.pyrograph_calcs as calc
import csv

def scatter_plot(quantity_id,data_set,quantities,groups,styles,output_directory,diagnostic_level):
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
    sigma_2_e = float(quantities[quantity_id]["Sigma_2_E"])/100.0
    
    x_data_set = []
    y_data_set = []
    
    # Margin size in inches.
    left_margin = 0.75
    right_margin = 0.25
    top_margin = 0.25
    bottom_margin = 0.75
    
    rc("figure.subplot", left=(left_margin/size))
    rc("figure.subplot", right=((size-right_margin)/size))
    rc("figure.subplot", bottom=(bottom_margin/size))
    rc("figure.subplot", top=((size-top_margin)/size))
    
    fig = plt.figure(figsize=(float(size),float(size))) # (w,h) values are size in inches.
    ax = fig.add_subplot(111, aspect='equal')
    
    # mu_sigma = calc.mu_2sigma(x_data_set,y_data_set,diagnostic_level)
    # mu_val = mu_sigma[0]
    # sigma_2_val = mu_sigma[1]
    # mu_max = plot_max*(1+mu_val)
    
    id_column = [quantity_id,title]
    
    tempfile = open('tempfile.csv', 'a')
    tempWriter = csv.writer(tempfile,delimiter=',')
    tempWriter.writerow(id_column)
    tempfile.close()
    
    # Draw Center Line
    ax.plot([plot_min,plot_max],[plot_min,plot_max], 'k-', linewidth=2.0, label='_nolegend_')
    
    # Iterate through data_set and plot by Group_Index.
    for group_id in data_set:
        group_title = groups[group_id]["Group_Title"]
        symboltype = styles[int(groups[group_id]["Style_Index"])]['Symbol_Style']
        symbolcolor = styles[int(groups[group_id]["Style_Index"])]['Fill_Color']
        symbolsize = int(styles[int(groups[group_id]["Style_Index"])]['Symbol_Size'])
        edgecolor = styles[int(groups[group_id]["Style_Index"])]['Edge_Color']
        
        #print "Data:",data_set[group_id]
        x_data = [x[0] for x in data_set[group_id]]
        y_data = [x[1] for x in data_set[group_id]]
        x_data_set.append(x_data)
        y_data_set.append(y_data)
        ax.scatter(x_data, y_data, s=symbolsize, c=symbolcolor, marker=symboltype, edgecolors=edgecolor, label=group_title)
    
    bias_sigma = calc.delta_sigma(x_data_set,y_data_set,sigma_2_e,diagnostic_level)
    bias = bias_sigma[0]
    sigma = bias_sigma[1]
    if diagnostic_level >= 3:
        print 'New Bias:',bias
        print 'New Sigma:',sigma
    bias_max = plot_max*bias
    
    line_style_list = [int(quantities[quantity_id]["Bias_Style"]),int(quantities[quantity_id]["Sigma_2_M_Style"]),int(quantities[quantity_id]["Sigma_2_E_Style"])]
    
    style_counter = 0
    pos = 0.
    space = 0.07

    for style_value in line_style_list:
        if style_value != 0:
            linecolor = styles[style_value]['Line_Color']
            linestyle = styles[style_value]['Line_Style']
            
            if styles[style_value]['Line_Width'] == 'None':
                print 'No Line Width specified in the style',style_value,'chosen for scatter plot statistic lines.'
                exit()
            else:
                linewidth = float(styles[style_value]['Line_Width'])
                
            if style_counter == 0:
                # Draw Bias Line
                bias_line = ax.plot([plot_min,plot_max],[plot_min,bias_max], c=linecolor, linestyle=linestyle, linewidth=linewidth, label='_nolegend_')
                pos = pos + space
                ax.text((plot_max-plot_min)*(title_position[0]+0.05), (plot_max-plot_min)*(title_position[1]-pos), r'$\mathrm{Bias}='+"%2.2f" % bias+'$')
                
            if style_counter == 1:
                # Upper 2 Sigma Line
                sigma2upper = ax.plot([plot_min,plot_max], [plot_min,(bias_max*(1+sigma*2))], c=linecolor, linestyle=linestyle, linewidth=linewidth, label='_nolegend_')
                # Lower 2 Sigma Line
                sigma2lower = ax.plot([plot_min,plot_max], [plot_min,(bias_max*(1-sigma*2))], c=linecolor, linestyle=linestyle, linewidth=linewidth, label='_nolegend_')
                pos = pos + space
                ax.text((plot_max-plot_min)*(title_position[0]+0.05), (plot_max-plot_min)*(title_position[1]-pos), r'$2 \, \widetilde{\sigma}_M='+"%2.2f" % (2*sigma)+'$')
                
            if style_counter == 2:
                # Draw Upper Exp Error
                exp_error_upper = ax.plot([plot_min,plot_max], [plot_min,(plot_max*(1+sigma_2_e))], c=linecolor, linestyle=linestyle, linewidth=linewidth, label='_nolegend_')
                # Draw Lower Exp Error
                exp_error_lower = ax.plot([plot_min,plot_max], [plot_min,(plot_max*(1-sigma_2_e))], c=linecolor, linestyle=linestyle, linewidth=linewidth, label='_nolegend_')
                pos = pos + space
                ax.text((plot_max-plot_min)*(title_position[0]+0.05), (plot_max-plot_min)*(title_position[1]-pos), r'$2 \, \widetilde{\sigma}_E='+"%2.2f" % (sigma_2_e)+'$')
                
            style_counter += 1
        else:
            style_counter += 1
            
    plt.xlabel(ind_axis_title)
    plt.ylabel(dep_axis_title)
    
    ax.text((plot_max-plot_min)*title_position[0], (plot_max-plot_min)*title_position[1], title, horizontalalignment='left') 
    
    ax.axis([plot_min, plot_max, plot_min, plot_max])
    
    if legend_loc != '':
        leg = ax.legend(loc=legend_loc)
        
    
    # if legend_loc != '':
    #                 for t in leg.get_texts():
    #                     t.set_fontsize('small')    # the legend text fontsize
            
    plt.show()
    plt.savefig(output_directory+output_filename)
    plt.close()

def comparison_plot(data_set,data_info,d1_index_set,d2_index_set,styles,output_directory,diagnostic_level):
    #Set Comparison Plot Settings
    title = data_info["Plot_Title"]
    size = float(data_info["Plot_Width"])
    output_filename = data_info["Plot_Filename"]
    legend_loc= data_info["Key_Position"]
    ind_axis_title = data_info["Ind_Title"]
    dep_axis_title = data_info["Dep_Title"]
    flip_axis = data_info["Flip_Axis"]
    xmin = float(data_info['Min_Ind'])
    xmax = float(data_info['Max_Ind'])
    ymin = float(data_info['Min_Dep'])
    ymax = float(data_info['Max_Dep'])
    title_position = eval(data_info["Title_Position"])
    title_x_pos = float(title_position[0])
    title_y_pos = float(title_position[1])
    
    mask_value = -9999.0
    
    # Margin size in inches.
    left_margin = 0.85
    right_margin = 0.15
    top_margin = 0.5
    bottom_margin = 0.5
    
    if data_info['d1_Key'][0] == '[':
        d1_column_names = eval(data_info['d1_Key'])
    else:
        d1_column_names = [data_info['d1_Key']]
        
    if data_info['d2_Key'][0] == '[':
        d2_column_names = eval(data_info['d2_Key'])
    else:
        d2_column_names = [data_info['d2_Key']]
    
    rc("figure.subplot", left=(left_margin/size))
    rc("figure.subplot", right=((size-right_margin)/size))
    rc("figure.subplot", bottom=(bottom_margin/(size*0.75)))
    rc("figure.subplot", top=(((size*0.75)-top_margin)/(size*0.75)))
    
    fig = plt.figure(figsize=(size,size*0.75)) # (w,h) values are size in inches.
    ax = fig.add_subplot(111)
    
    plot_type_function = {'linear':ax.plot,'loglog':ax.loglog,'semilogx':ax.semilogx,'semilogy':ax.semilogy}
    
    d1_data = data_set[0]
    d2_data = data_set[1]
    d1_ind_data = d1_data[0]
    d2_ind_data = d2_data[0]
    
    d1_dep_data = d1_data[1:]
    d2_dep_data = d2_data[1:]
    
    # Mask out data set to -9999.0 in y_data_set
    d1_dep_data_to_plot = []
    d2_dep_data_to_plot = []
    
    for d1_dep_data_subset in d1_dep_data:
        dep_values = M.array(d1_dep_data_subset)
        d1_dep_values_masked = M.masked_where(dep_values == mask_value , dep_values)
        d1_dep_data_to_plot.append(d1_dep_values_masked)
    
    for d2_dep_data_subset in d2_dep_data:
        dep_values = M.array(d2_dep_data_subset)
        d2_dep_values_masked = M.masked_where(dep_values == mask_value , dep_values)
        d2_dep_data_to_plot.append(d2_dep_values_masked)
    
    linecount = 0
    for d1_dep_masked_data in d1_dep_data_to_plot:
        
        d1_style = eval(data_info["d1_Style"])
        d1_key = d1_column_names[linecount]
        symboltype = styles[int(d1_style[linecount])]['Symbol_Style']
        edgecolor = styles[int(d1_style[linecount])]['Edge_Color']
        symbolcolor = styles[int(d1_style[linecount])]['Fill_Color']
        
        if styles[int(d1_style[linecount])]['Symbol_Size'] == 'None':
            symbolsize = 'None'
        else:
            symbolsize = int(styles[int(d1_style[linecount])]['Symbol_Size'])
        
        linecolor = styles[int(d1_style[linecount])]['Line_Color']
        linestyle = styles[int(d1_style[linecount])]['Line_Style']
        
        if styles[int(d1_style[linecount])]['Line_Width'] == 'None':
            linewidth = 'None'
        else:
            linewidth = float(styles[int(d1_style[linecount])]['Line_Width'])
            
        if flip_axis == 'yes':
            d1_x_data = d1_dep_masked_data[d1_index_set[0]:d1_index_set[1]+1]
            d1_y_data = d1_ind_data[d1_index_set[0]:d1_index_set[1]+1]
        else:
            d1_x_data = d1_ind_data[d1_index_set[0]:d1_index_set[1]+1]
            d1_y_data = d1_dep_masked_data[d1_index_set[0]:d1_index_set[1]+1]
        
        plot_type_function[data_info['Plot_Type']](d1_x_data, d1_y_data, c=linecolor, linestyle=linestyle, linewidth=linewidth, marker=symboltype, ms=symbolsize, mfc=symbolcolor, mec=edgecolor, label=d1_key)
        linecount += 1
        
    linecount = 0
    for d2_dep_masked_data in d2_dep_data_to_plot:
        
        d2_style = eval(data_info["d2_Style"])
        d2_key = d2_column_names[linecount]
        symboltype = styles[int(d2_style[linecount])]['Symbol_Style']
        edgecolor = styles[int(d2_style[linecount])]['Edge_Color']
        symbolcolor = styles[int(d2_style[linecount])]['Fill_Color']
        
        if styles[int(d2_style[linecount])]['Symbol_Size'] == 'None':
            symbolsize = 'None'
        else:
            symbolsize = int(styles[int(d2_style[linecount])]['Symbol_Size'])
        
        linecolor = styles[int(d2_style[linecount])]['Line_Color']
        linestyle = styles[int(d2_style[linecount])]['Line_Style']
        
        if styles[int(d2_style[linecount])]['Line_Width'] == 'None':
            linewidth = 'None'
        else:
            linewidth = float(styles[int(d2_style[linecount])]['Line_Width'])
            
        if flip_axis == 'yes':
            d2_x_data = d2_dep_masked_data[d2_index_set[0]:d2_index_set[1]+1]
            d2_y_data = d2_ind_data[d2_index_set[0]:d2_index_set[1]+1]
        else:
            d2_x_data = d2_ind_data[d2_index_set[0]:d2_index_set[1]+1]
            d2_y_data = d2_dep_masked_data[d2_index_set[0]:d2_index_set[1]+1]
        
        plot_type_function[data_info['Plot_Type']](d2_x_data, d2_y_data, c=linecolor, linestyle=linestyle, linewidth=linewidth, marker=symboltype, ms=symbolsize, mfc=symbolcolor, mec=edgecolor, label=d2_key)
        linecount += 1
    
    if flip_axis == 'yes':
        plt.xlabel(dep_axis_title)
        plt.ylabel(ind_axis_title)
        
    else:
        plt.xlabel(ind_axis_title)
        plt.ylabel(dep_axis_title)
        
    if flip_axis == 'yes':
        if data_info['Plot_Type'] == 'linear':
            ax.text(ymin+(ymax-ymin)*title_x_pos, xmin+(xmax-xmin)*title_y_pos, title, horizontalalignment='left')
        elif data_info['Plot_Type'] == 'loglog':
            pos1 = 10**(log10(ymin)+((log10(ymax)-log10(ymin))*title_x_pos))
            pos2 = 10**(log10(xmin)+((log10(xmax)-log10(xmin))*title_y_pos))
            ax.text(pos1, pos2, title, horizontalalignment='left')
        elif data_info['Plot_Type'] == 'semilogx':
            pos1 = 10**(log10(ymin)+((log10(ymax)-log10(ymin))*title_x_pos))
            pos2 = xmin+(xmax-xmin)*title_y_pos
            ax.text(pos1, pos2, title, horizontalalignment='left')
        elif data_info['Plot_Type'] == 'semilogy':
            pos1 = ymin+(ymax-ymin)*title_x_pos
            pos2 = 10**(log10(xmin)+((log10(xmax)-log10(xmin))*title_y_pos))
            ax.text(pos1, pos2, title, horizontalalignment='left')
    else:
        if data_info['Plot_Type'] == 'linear':
            ax.text(xmin+(xmax-xmin)*title_x_pos, ymin+(ymax-ymin)*title_y_pos, title, horizontalalignment='left')
        elif data_info['Plot_Type'] == 'loglog':
            pos1 = 10**(log10(xmin)+((log10(xmax)-log10(xmin))*title_x_pos))
            pos2 = 10**(log10(ymin)+((log10(ymax)-log10(ymin))*title_y_pos))
            ax.text(pos1, pos2, title, horizontalalignment='left')
        elif data_info['Plot_Type'] == 'semilogx':
            pos1 = 10**(log10(xmin)+((log10(xmax)-log10(xmin))*title_x_pos))
            pos2 = ymin+(ymax-ymin)*title_y_pos
            ax.text(pos1, pos2, title, horizontalalignment='left')
        elif data_info['Plot_Type'] == 'semilogy':
            pos1 = xmin+(xmax-xmin)*title_x_pos
            pos2 = 10**(log10(ymin)+((log10(ymax)-log10(ymin))*title_y_pos))
            ax.text(pos1, pos2, title, horizontalalignment='left')
    
    if flip_axis == 'yes':
        ax.axis([ymin,ymax,xmin,xmax])
    else:
        ax.axis([xmin,xmax,ymin,ymax])
    
    if legend_loc != 'none':
        leg = ax.legend(loc=legend_loc)
    
    # if legend_loc != 'none':
    #                     for t in leg.get_texts():
    #                         t.set_fontsize('small')    # the legend text fontsize
            
    plt.show()
    plt.savefig(output_directory+output_filename)
    plt.close()


# def parallel_plot(plot_type,plot_args,num_cpus=1):
#     import math, sys, time
#     import pp
#     
#     # Tuple of all parallel python servers to connect with
#     ppservers = ()
#     #ppservers = ("10.0.0.1",)
#     
#     if num_cpus >= 1:
#         # Creates jobserver with ncpus workers
#         job_server = pp.Server(num_cpus, ppservers=ppservers)
#     else:
#         # Creates jobserver with automatically detected number of workers
#         job_server = pp.Server(ppservers=ppservers)
# 
#     print "Starting parallel plotting with", job_server.get_ncpus(), "workers."
#     
#     if plot_type == 's':
#         #print "Running Scatter Plots"
#         func = scatter_plot()
#     elif plot_type == 'c':
#         #print "Running Comparison Plots"
#         func = comparison_plot()
#         
#     else:
#         print "Plot type is not recognized."
#         exit()
#     
#     start_time = time.time()
    
    
