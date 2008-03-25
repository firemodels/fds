import csv
from re import *
from pyx import *


### Set Global Variables

data_directory = "../../Validation/"
output_directory = "../../Manuals/FDS_5_Validation_Guide/FIGURES/"
config_file_name = "Validation_Data_Config_File.csv"

scatter_data_dict = {}
combined_scatter_data = {}

### Define Functions

def extract_config_data(config_file):
    # Collect Data from Config File
    quantities_dict = {}
    data_dict = {}
    keyed_quantities = {}
    keyed_data = {}
    quantity_counter = 0
    data_counter = 0
    
    fh = file(config_file, 'U')
    data_array = csv.reader(fh)
    
    config_lists = [list(sublist) for sublist in data_array]
    #print config_lists
    print str(len(config_lists)) + " lines read from Configuration file\n"
    
    #Build Quantity and Data Dictionaries, with values keyed by config column name.
    for list_item in config_lists:
        if list_item[0] == 'q':
            if quantity_counter < 1:
                quantity_header = list_item[2:]
            if quantity_counter >= 1:
                for x in range(len(quantity_header)):
                    keyed_quantities[quantity_header[x]] = list_item[x+2]
                quantities_dict[list_item[1]] = keyed_quantities
                keyed_quantities = {}
            quantity_counter += 1
        elif list_item[0] == 'd':
            if data_counter < quantity_counter:
                data_counter = quantity_counter - 1
                data_header = list_item[1:]
            if data_counter >= quantity_counter:
                for x in range(len(data_header)):
                    keyed_data[data_header[x]] = list_item[x+1]
                data_key_name = list_item[2].strip()+"-"+list_item[6].strip()
                print data_key_name
                data_dict[data_key_name] = keyed_data
                keyed_data = {}
            data_counter += 1
        else:
            print """No d or q, skip row."""
    
    #Convert quantities keys from string to integer.
    tempdict = dict((int(key), values) for key,values in quantities_dict.iteritems())
    quantities_dict = tempdict
    
    # Return a single list object containing the dictionaries.
    return [quantities_dict,data_dict]

def compute_difference(m_p,m_o,e_p,e_o):
    # Equation to compute relative difference of model and experimental data.
    M_p = float(m_p) #Peak value of Model Prediction to float
    M_o = float(m_o) #Original value of Model Prediction to float
    E_p = float(e_p) #Peak value of Experimental Measurement to float
    E_o = float(e_o) #Original value of Experimental Measurement to float
    e = (((M_p-M_o)-(E_p-E_o))/(E_p-E_o))
    return e

def extract_comp_data(comp_file_info):
    ## Read in from config file and Process data from Source .csv files.
    
    exp_data_dict = {}
    mod_data_dict = {}
    
    #List of variables from configuration file column names.
    
    exp_data_filename = comp_file_info['Exp_Filename'] #String of filename
    exp_column_name = comp_file_info['Exp_Col_Name'].strip() #Experimental Data Column Name
    exp_column_name_row_num = int(comp_file_info['Exp_Col_Name_Row'])-1 #Experimental Data Column Name Row Number
    exp_data_row_num = int(comp_file_info['Exp_Data_Row'])-1 #Experimental Data Starting Row Number
    
    mod_data_filename = comp_file_info['Mod_Filename'] #String of filename
    mod_column_name = comp_file_info['Mod_Col_Name'].strip() #Modeling Data Column Name
    mod_column_name_row_num = int(comp_file_info['Mod_Col_Name_Row'])-1 #Modeling Data Column Name Row Number
    mod_data_row_num = int(comp_file_info['Mod_Data_Row'])-1 #Modeling Data Starting Row Number
    
    scatter_data_label = comp_file_info['Quantity']+"-"+comp_file_info['Dataname']+"-"+comp_file_info['Exp_Col_Name']
    
    exp_start_time_data_val = comp_file_info['Exp_Start_(min.)'] #String in minutes to start exp plot data
    exp_stop_time_data_val = comp_file_info['Exp_End_(min.)']  #String in minutes to stop exp plot data
    exp_start_time_comp_val = comp_file_info['Exp_Comp_Start_(min.)'] #String in minutes to start exp compare data
    exp_stop_time_comp_val = comp_file_info['Exp_Comp_End_(min.)']  #String in minutes to start exp compare data
    exp_initial_value = comp_file_info['Exp_Intitial_Value']       #Initial Value for Quantity
    
    mod_start_time_data_val = comp_file_info['Mod_Start_(min.)'] #String in minutes to start mod plot data
    mod_stop_time_data_val = comp_file_info['Mod_End_(min.)']  #String in minutes to stop mod plot data
    mod_start_time_comp_val = comp_file_info['Mod_Comp_Start_(min.)'] #String in minutes to start mod compare data
    mod_stop_time_comp_val = comp_file_info['Mod_Comp_End_(min.)']  #String in minutes to start mod compare data
    mod_initial_value = comp_file_info['Mod_Intitial_Value']       #Initial Value for Quantity
    
    min_max = comp_file_info['max/min'] #String indicating if min or max value is required.

    def find_start_stop_index(data_dict,col_name,start_time_data,stop_time_data,start_time_comp,stop_time_comp):
        #This function is used to find index numbers for start and stop points in plotting and min-max values.
        
        rowcounter1 = 0
        for time_value1 in data_dict[col_name]:
            if time_value1 >= (float(start_time_data)*60):
                #print "Set #1"
                #print "Time Starts at row #:", str(rowcounter1)
                #print "With a value of:", str(time_value1)
                time_start_index = rowcounter1
                break                
            rowcounter1 += 1
        
        rowcounter2 = 0
        for time_value2 in data_dict[col_name]:
            if float(data_dict[col_name][(len(data_dict[col_name])-1)]) < (float(stop_time_data)*60):
                #print "Specified end of plot time is greater than end of time in the data set. \nUsing last value in the time column.\n"
                #print "Time used is: "+str(float(data_dict[col_name][(len(data_dict[col_name])-1)]))+"\n"
                time_end_index = (len(data_dict[col_name])-1)
                break
            else:
                row_number2 = (rowcounter2 - 1)
                #print "Set #2"
                #print "Time Ends at row #: "+str(row_number2)
                #print "With a value of: "+str(data_dict[col_name][row_number2])
                time_end_index = row_number2
                break

            if time_value2 < (float(stop_time_data)*60):
                rowcounter2 += 1
                
        rowcounter3 = 0                
        for time_value3 in data_dict[col_name]:
            if time_value3 >= (float(start_time_comp)*60):
                #print "Set #3"
                #print "Comparison Time Starts at row #:", str(rowcounter3)
                #print "With a value of:", str(time_value3)
                minmax_start_index = rowcounter3
                break
            rowcounter3 += 1
        
        rowcounter4 = 0
        for time_value4 in data_dict[col_name]:
            if float(data_dict[col_name][(len(data_dict[col_name])-1)]) < (float(stop_time_comp)*60):
                #print "Specified end of comparison time is greater than end of time in the data set. \nUsing last value in the time column."
                #print "Time used is: "+str(float(data_dict[col_name][(len(data_dict[col_name])-1)]))+"\n"
                minmax_end_index = (len(data_dict[col_name])-1)
                break

            if time_value4 < (float(stop_time_data)*60):
                rowcounter4 += 1
            else:
                row_number4 = (rowcounter4 - 1)
                #print "Set #4"
                #print "Comparison Time Ends at row #: "+str(row_number4)
                #print "With a value of: "+str(data_dict[col_name][row_number4])
                minmax_end_index = row_number4
                break

        return (time_start_index, time_end_index, minmax_start_index, minmax_end_index)

    exp_file_object = open(data_directory+exp_data_filename, "U")
    mod_file_object = open(data_directory+mod_data_filename, "U")
    
    ## Start File Processing
    
    #Read in experimental data and flip lists from rows to columns.
    print "Reading in:", exp_data_filename
    exp_data_cols = zip(*csv.reader(exp_file_object))
    #Convert tuples to lists.
    exp_data_list = [list(sublist) for sublist in exp_data_cols]
    #Pull the Time column name out and strip whitespace from ends of string.
    exp_time_col_name = exp_data_list[0][exp_column_name_row_num].strip()
    
    #Build Experimental Data Dictionary
    for exp_list in exp_data_list:
        exp_data_dict[exp_list[exp_column_name_row_num].strip()] = map(float, exp_list[exp_data_row_num:])
        #print "Exp. Data Dict:", exp_data_dict[(exp_list[exp_column_name_row_num]).strip()]
    
    #Read in model data and flip lists from rows to columns.
    print "Reading in:", mod_data_filename
    mod_data_cols = zip(*csv.reader(mod_file_object))
    #Convert tuples to lists.
    mod_data_list = [list(sublist) for sublist in mod_data_cols]
    #Pull the Time column name out and strip whitespace from ends of string.
    mod_time_col_name = mod_data_list[0][mod_column_name_row_num].strip()
    
    #Build Prediction/Model Data Dictionary
    for mod_list in mod_data_list:
        mod_data_dict[mod_list[mod_column_name_row_num].strip()] = map(float, mod_list[mod_data_row_num:])
        #print "Model Data Dict:", mod_data_dict[(mod_list[mod_column_name_row_num]).strip()]

    exp_comp_ranges = find_start_stop_index(exp_data_dict,exp_time_col_name,exp_start_time_data_val,exp_stop_time_data_val,exp_start_time_comp_val,exp_stop_time_comp_val)
    #print exp_comp_ranges
    
    mod_comp_ranges = find_start_stop_index(mod_data_dict,mod_time_col_name,mod_start_time_data_val,mod_stop_time_data_val,mod_start_time_comp_val,mod_stop_time_comp_val)
    #print mod_comp_ranges
    
    ##Find max or min values.
    
    #Exp min_max value
    exp_data_values_comp = exp_data_dict[exp_column_name][exp_comp_ranges[2]:exp_comp_ranges[3]]

    if  min_max == 'max':
        #print min_max, str(max(exp_data_values_comp))
        exp_peak_value = max(exp_data_values_comp)
    elif min_max == 'min':
        #print min_max, str(min(exp_data_values_comp))
        exp_peak_value = min(exp_data_values_comp)
    else:
        print "Min or Max is undefined in the input file."
    
    #Mod min_max value
    mod_data_values_comp = mod_data_dict[mod_column_name][mod_comp_ranges[2]:mod_comp_ranges[3]]
    #print mod_data_values_comp
    
    if  min_max == 'max':
        #print min_max, str(max(mod_data_values_comp))
        mod_peak_value = max(mod_data_values_comp)
    elif min_max == 'min':
        #print min_max, str(min(mod_data_values_comp))
        mod_peak_value = min(mod_data_values_comp)
    else:
        print "Min or Max is undefined in the input file."
    
    #print mod_peak_value
    #print mod_initial_value
    #print exp_peak_value
    #print exp_initial_value
    
    # This allows the d line Quantity value to be set to 0 when either model or experimental data is missing.
    if comp_file_info['Quantity'] == 0:
        print "Quantity set to 0, no comparison made."
        relative_difference = 'NA'
    else:
        relative_difference = compute_difference(mod_peak_value,mod_initial_value,exp_peak_value,exp_initial_value)
    
    #Append Min_Max Values to Global Scatter Data Dictionary.
    scatter_data_dict[scatter_data_label] = [exp_peak_value,mod_peak_value,relative_difference]
    
    #Create data lists based on specified ranges
    exp_data_seconds = zip(exp_data_dict[exp_time_col_name][exp_comp_ranges[0]:exp_comp_ranges[1]], exp_data_dict[exp_column_name][exp_comp_ranges[0]:exp_comp_ranges[1]])
    #print exp_data_seconds
    
    mod_data_seconds = zip(mod_data_dict[mod_time_col_name][mod_comp_ranges[0]:mod_comp_ranges[1]], mod_data_dict[mod_column_name][mod_comp_ranges[0]:mod_comp_ranges[1]])
    #print mod_data_seconds
    
    #Convert time to minutes from seconds.
    exp_data = [[x[0] / 60, x[1]] for x in exp_data_seconds]
    #print exp_data
    mod_data = [[x[0] / 60, x[1]] for x in mod_data_seconds]
    #print mod_data
    
    
    # Return list of X,Y lists.
    #exp_data=[[exp_time,exp_quantity_value],[exp_time,exp_quantity_value]]
    #mod_data=[[mod_time,mod_quantity_value],[mod_time,mod_quantity_value]]
    
    return [exp_data,mod_data]


## Plot Validation Data to PDF files
## Two kinds of plots... Comparison and Scatter

def comparison_plot(plot_data,exp_data,mod_data):
    #plot_data is the items from the 'd' rows of the config file.
    #Variables from configuration file column names.

    
    # Variables for plot.
    plot_title = plot_data['Plot_Title']
    #print plot_title
    x_title = plot_data['X_Title']
    y_title = plot_data['Y_Title']   
    min_x = float(plot_data['Min_X'])
    max_x = float(plot_data['Max_X'])
    min_y = float(plot_data['Min_Y'])
    max_y = float(plot_data['Max_Y'])
    title_quadrant = int(plot_data['Title_Quadrant'])
    key_pos = plot_data['Key_Position']
    key_dist = 0.2*unit.v_cm
    plot_width = int(plot_data['Plot_Width(cm)'])
    
    #Set plot legend key text.
    exp_key = plot_data['Exp_Key']
    mod_key = plot_data['Mod_Key']
    
    #Create filename from fields in input file record.
    plot_file_name = plot_data["Plot_Filename"]
    
    # Determine the location for the key, alignment based on key_quadrant setting.
    # Replace quad code with actual position letters
    if key_pos == "tl" or "tc" or "tr" or "ml" or "mc" or "mr" or "bl" or "bc" or "br":
        ()
        #print "Key Position =", key_pos
    else:
        print "The key position was not specified./nUsing the default bottom right position."
        key_pos = "br"
    
    #Begin Plotting
    # Initialize graph object
    g = graph.graphxy(width=plot_width, ratio=4./3, key=graph.key.key(pos=key_pos, dist=key_dist), 
                        x=graph.axis.linear(title=x_title, min=min_x, max=max_x), 
                        y=graph.axis.linear(title=y_title, min=min_y, max=max_y))

    #One line at a time added to plot from each data set.
    # Plot Experimental data
    g.plot(graph.data.points(exp_data, title=exp_key, x=1, y=2),
        [graph.style.line([style.linewidth.Thick, style.linestyle.solid])])
    # Plot Predicted/Model data
    g.plot(graph.data.points(mod_data, title=mod_key, x=1, y=2),
        [graph.style.line([color.rgb.black, style.linewidth(0.06*unit.w_cm), style.linestyle.dotted])])
    
    # Now plot the Title text, alignment based on title quadrant setting.
    if title_quadrant == 1:
        g.text(0.1, g.height - 0.2, plot_title, [text.halign.left, text.valign.top, text.size.small])
    elif title_quadrant == 2:
        g.text(g.width-0.1, g.height - 0.2, plot_title, [text.halign.right, text.valign.top, text.size.normalsize])
    elif title_quadrant == 3:
        g.text(0.1, 0.2, plot_title, [text.halign.left, text.valign.bottom, text.size.normalsize])
    elif title_quadrant == 4:
        g.text(g.width-0.1, 0.2, plot_title, [text.halign.right, text.valign.bottom, text.size.normalsize])
    else:
        print "A quadrant for the title location was not specified./nUsing the default top left quadrant."
        g.text(0.1, g.height - 0.2, plot_title, [text.halign.left, text.valign.top, text.size.small])
        
    # Write the output
    plot_file_path = output_directory+plot_file_name
    g.writePDFfile(plot_file_path)
    print "Plot to: \n", plot_file_path+".PDF"


def scatter_plot(plot_info,data_set):
    #plot_info is details about the overall plot layout and text.
    #data_set is a dictionary keyed by test name containing lists of X and Y data points.
    print plot_info
    #print data_set
    
    # Need quantity and group iterations.  Also need way to make key based on groups.
    for quantity_number in plot_info:
        if data_set[quantity_number] == []:
            print "No Scatter Plot Data in Quantity "+str(quantity_number)+" Dataset."
        else:
            # Set variables for Plot extracted from the first group of lines in config file starting with 'q'.
            
            # Variables for plot.
            plot_title = plot_info[int(quantity_number)]['Scatter_Plot_Title']
            print plot_title
            x_title = plot_info[int(quantity_number)]['X_Title']
            y_title = plot_info[int(quantity_number)]['Y_Title']
            min_x = float(plot_info[int(quantity_number)]['Plot_Min'])
            #print min_x
            max_x = float(plot_info[int(quantity_number)]['Plot_Max'])
            #print max_x
            min_y = float(plot_info[int(quantity_number)]['Plot_Min'])
            max_y = float(plot_info[int(quantity_number)]['Plot_Max'])
            percent_error = int(plot_info[int(quantity_number)]['%error'])
            title_quadrant = int(plot_info[int(quantity_number)]['Title_Quadrant'])
            key_pos = plot_info[int(quantity_number)]['Key_Position']
            key_dist = 0.2*unit.v_cm
            plot_width = int(plot_info[int(quantity_number)]['Plot_Width(cm)'])
            
            #Create filename from fields in input file record.
            plot_file_name = plot_info[int(quantity_number)]['Plot_Filename']
            #print plot_file_name
            
            # Determine the location for the key, alignment based on key_quadrant setting.
            if key_pos == "tl" or "tc" or "tr" or "ml" or "mc" or "mr" or "bl" or "bc" or "br":
                ()
                #print "Key Position =", key_pos
            else:
                print "The key position was not specified./nUsing the default bottom right position."
                key_pos = "br"
            
            #Begin Plotting
            #print exp_data
            #print mod_data
            # Initialize graph object
            g = graph.graphxy(width=plot_width, ratio=1/1, key=graph.key.key(pos=key_pos, dist=key_dist), 
                                x=graph.axis.linear(title=x_title, min=min_x, max=max_x), 
                                y=graph.axis.linear(title=y_title, min=min_y, max=max_y))
            
            #Plot Midline and Error bounds lines.
            errorLineCenterPoints = [[min_x,min_y],[max_x,max_y]]
            #print errorLineCenterPoints
            lower_bound = max_y - max_y * percent_error / 100
            #print lower_bound
            errorLineLowerPoints = [[min_x,min_y],[max_x,lower_bound]]
            #print errorLineLowerPoints
            upper_bound = max_y + max_y * percent_error / 100.0
            #print upper_bound
            errorLineUpperPoints = [[min_x,min_y],[max_x,upper_bound]]
            #print errorLineUpperPoints
            
            g.plot(graph.data.points(errorLineCenterPoints, title=None, x=1, y=2),
                    [graph.style.line([style.linewidth.Thin, style.linestyle.solid])])
                    
            g.plot(graph.data.points(errorLineLowerPoints, title=None, x=1, y=2),
                    [graph.style.line([style.linewidth.Thin, style.linestyle.dashed])])
                    
            g.plot(graph.data.points(errorLineUpperPoints, title=None, x=1, y=2),
                    [graph.style.line([style.linewidth.Thin, style.linestyle.dashed])])
        
            #One line at a time added to plot from each data set.
            # Iterate over items in data dictionary key for keys that are not [].
            # Append data sets to scatter_plot_data_list
            scatter_plot_data_list = []
            mystyle = graph.style.symbol(
                            graph.style.symbol.changetriangletwice, 
                            size=0.1*unit.v_cm, 
                            symbolattrs=[graph.style.symbol.changefilledstroked, 
                                        attr.changelist([color.rgb.red, color.rgb.green, color.rgb.blue])])
            for data in data_set[quantity_number]:
                #IF data[1] is a list of lists then loop over each group of data with a different symbol.
                #ELSE 
                scatter_plot_data = [data[1]]
                #print scatter_plot_data
                g.plot(graph.data.points(scatter_plot_data, x=1, y=2, title=data[0].replace('_',' ')), [mystyle])

            # Now plot the Title text, alignment based on title quadrant setting.
            if title_quadrant == 1:
                g.text(0.1, g.height - 0.2, plot_title, [text.halign.left, text.valign.top, text.size.small])
            elif title_quadrant == 2:
                g.text(g.width-0.1, g.height - 0.2, plot_title, [text.halign.right, text.valign.top, text.size.normalsize])
            elif title_quadrant == 3:
                g.text(0.1, 0.2, plot_title, [text.halign.left, text.valign.bottom, text.size.normalsize])
            elif title_quadrant == 4:
                g.text(g.width-0.1, 0.2, plot_title, [text.halign.right, text.valign.bottom, text.size.normalsize])
            else:
                print "A quadrant for the title location was not specified./nUsing the default top left quadrant."
                g.text(0.1, g.height - 0.2, plot_title, [text.halign.left, text.valign.top, text.size.small])
                
            # Write the output
            plot_file_path = output_directory+plot_file_name
            #print plot_file_path
            g.writePDFfile(plot_file_path)
            print "Plot to: \n", plot_file_path+".PDF"
    

### Start of Main Code

#Get information from config file.
config_and_data_dicts = extract_config_data(config_file_name)
print "There are "+str(len(config_and_data_dicts[0]))+" quantities to process."
print "There are "+str(len(config_and_data_dicts[1]))+" data sets to plot."

## Create comparison plots
for data_record in config_and_data_dicts[1]:
    # Extract relevant portions of comparison data as defined in config file.
    comp_data_to_plot = extract_comp_data(config_and_data_dicts[1][data_record])
    
    #Seperate experimental and model data lists.
    exp_plot_data = comp_data_to_plot[0]
    mod_plot_data = comp_data_to_plot[1]
    
    # Create plot for data_record.
    comparison_plot(config_and_data_dicts[1][data_record],exp_plot_data,mod_plot_data)


## Create scatter plots
scatter_quantity = 1
temp_scatter_data_list = []

#Grouping Scatter Data by Quantity
for scatter_plot_record in sorted(config_and_data_dicts[0]):
    #print "Scatter record:", scatter_plot_record
    #print "Quantity:", config_and_data_dicts[0][scatter_plot_record]['Comparison_Quantities']
    for scatter_data_key in sorted(scatter_data_dict):
        split_key = split("-",scatter_data_key)
        if scatter_data_key[0] == str(scatter_plot_record) and scatter_plot_record != []:
            #print "Test name:", split_key[1]
            #print "Scatter Data:", scatter_data_dict[scatter_data_key][:2]
            temp_scatter_data_list.append([split_key[1], scatter_data_dict[scatter_data_key][:2]])
            
    combined_scatter_data[scatter_quantity] = temp_scatter_data_list
    temp_scatter_data_list = []
    scatter_quantity = scatter_quantity + 1
        
    #scatter_plot_info = config_and_data_dicts[1][scatter_plot_record]
    #print scatter_plot_info    print scatter_data_key, scatter_data_dict[scatter_data_key]

# Create plots
scatter_plot(config_and_data_dicts[0],combined_scatter_data)


## Write Summary Data to File.
#NRC Comparisons Output

# Output for each data set 
#*Exp Zero Val
#*Exp Peak Val
#*Peak Time Val
#*Mod Zero Val
#*Mod Peak Val
#*Peak Time Val
#*DeltaE
#*DeltaM
#*Rel Diff
