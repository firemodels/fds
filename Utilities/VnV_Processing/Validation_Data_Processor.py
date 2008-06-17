import csv
from re import *
from pyx import *

### Gather Processing Preferences from User
## Validation, Verification or BOTH
## Plot Comparison, Scatter or BOTH

### Set Global Variables for Validation

data_directory = "../../Validation/"
output_directory = "../../Manuals/FDS_5_Validation_Guide/FIGURES/"
config_file_name = "Validation_Data_Config_File.csv"

### Set Global Variables for Verification

#data_directory = "../../Verification/"
#output_directory = "../../Manuals/FDS_5_User_Guide/FIGURES/"
#config_file_name = "Verification_Data_Config_File.csv"

### Set General Global Variables
scatter_data_dict = {}
combined_scatter_data = {}

### Define Functions

def extract_config_data(config_file):
    print "*** Extracting Configuration Data and Building Dictionaries ***"
    
    # Collect Data from Config File
    groups_dict = {}
    quantities_dict = {}
    data_dict = {}
    keyed_groups = {}
    keyed_quantities = {}
    keyed_data = {}
    quantity_counter = 0
    data_counter = 0
    group_counter = 0
    skip_counter = 0
    
    #Create File Object
    try:
        fh = file(config_file, 'U')
    except:
	    print"!!! The Config File "+config_file_name+" does not exist or the path defined in the script is incorrect. !!!"
	    exit()
	
    #Read file with csv module.
    data_array = csv.reader(fh)
    #Convert into List object
    config_lists = [list(sublist) for sublist in data_array]
    print str(len(config_lists))+" lines read in from "+config_file_name+"\n"
    
    #Build Quantity and Data Dictionaries, with values keyed by config column name.
    for list_item in config_lists:
        if list_item[0] == 'g':
            if group_counter < 1:
                group_header = list_item[2:]
            else:
                if group_counter >= 1:
                    for x in range(len(group_header)):
                        keyed_groups[group_header[x]] = list_item[x+2]
                    groups_dict[int(list_item[1])] = keyed_groups
                    keyed_groups = {}
            group_counter += 1
            #print "This is group line #"+str(group_counter)+"."
        elif list_item[0] == 'q':
            if quantity_counter < 1:
                quantity_header = list_item[2:]
            if quantity_counter >= 1:
                for x in range(len(quantity_header)):
                    keyed_quantities[quantity_header[x]] = list_item[x+2]
                #print "List item 1:",int(list_item[1])
                quantities_dict[int(list_item[1])] = keyed_quantities
                keyed_quantities = {}
            quantity_counter += 1
        elif list_item[0] == 'd':
            if data_counter < quantity_counter:
                data_counter = quantity_counter - 1
                data_header = list_item[1:]
            if data_counter >= quantity_counter:
                for x in range(len(data_header)):
                    keyed_data[data_header[x].strip()] = list_item[x+1]
                data_key_name = keyed_data['Quantity'].strip()+"~"+keyed_data['Group'].strip()+"~"+keyed_data['Dataname'].strip()+"~"+keyed_data['Exp_Y_Col_Name'].strip()
                #print "Key Name:", data_key_name
                data_dict[data_key_name] = keyed_data
                #print data_dict[data_key_name]
                keyed_data = {}
            data_counter += 1
        else:
            skip_counter = skip_counter + 1
            #print """No g, d or q, skip row."""
            
    print "There were "+str(skip_counter)+" lines skipped, out of the "+str(len(config_lists))+" lines read in."
    # Return a single list object containing the dictionaries.
    #print groups_dict
    #print quantities_dict
    #print data_dict
    return [groups_dict,quantities_dict,data_dict]

def find_start_stop_index(data_dict,col_name,start_data,stop_data,start_comp,stop_comp,x_scale):
    #This function is used to find index numbers for start and stop points in plotting and min-max values.
    print "\n*** Find Start/End Indexes for X Column ***"
    #print "X Scale Factor:", x_scale
    #print "Stop Data:", stop_data
    print "X Column name:", col_name
    #print "X Column Data:", data_dict[col_name]
    rowcounter1 = 0
    print "*** Finding X Column Start Index ***"
    for value1 in data_dict[col_name]:
        #print "Value 1:", value1, " : ", rowcounter1
        if value1 == 'Null':
            continue
        else:
            if value1 >= (float(start_data)*float(x_scale)):
                #print "X Column Starts at row #:", str(rowcounter1)
                #print "With a value of:", str(value1)
                Xcol_start_index = rowcounter1
                break
        rowcounter1 = rowcounter1 + 1
    rowcounter2 = 0
    end_index = 0
    index_count = 0
    print "*** Finding Index of last good value in X Column ***"    
    for end_data_index in data_dict[col_name]:
        #print end_data_index, " : ", index_count
        if end_data_index == 'Null':
            end_index = index_count - 1
            break
        else:
            index_count = index_count + 1
    if end_index == 0:
        end_index = index_count - 1
    #print "Length of X Column:", end_index + 1
    #print "Last value in X Column:", data_dict[col_name][end_index]
    #print "Index of Last Numeric Value:", end_index
    
    print "*** Finding X Column End Index ***"
    for value2 in data_dict[col_name]:
        #print "Value 2:", value2, " : ", rowcounter2
        if float(data_dict[col_name][end_index]) < (float(stop_data)*float(x_scale)):
            print "Specified end of X Column Data is greater than end of time in the data set. \nUsing last value in the time column.\n"
            print "Value used is: "+str(float(data_dict[col_name][end_index]))+"\n"
            Xcol_end_index = (end_index)
            break
        if value2 < (float(stop_data)*float(x_scale)):
            rowcounter2 = rowcounter2 + 1
        else:
            row_number2 = (rowcounter2)
            print "X Column Ends at Index #: ", str(row_number2), "with a value of: ", str(data_dict[col_name][row_number2])
            Xcol_end_index = row_number2
            break
    rowcounter3 = 0
    print "*** Finding X Column Comp Start Index ***"
    for value3 in data_dict[col_name]:
        if value3 >= (float(start_comp)*float(x_scale)):
            print "Comparison Starts at Index #:", str(rowcounter3), "with a value of:", str(value3)
            minmax_start_index = rowcounter3
            break
        rowcounter3 = rowcounter3 + 1 
    
    rowcounter4 = 0
    print "*** Finding X Column Comp End Index ***"
    for value4 in data_dict[col_name]:
        scaled_stop_value = (float(stop_comp)*float(x_scale))
        end_index_value = float(data_dict[col_name][end_index])
        #print "Scaled Stop Value and End Index Value:", scaled_stop_value, " : ",end_index_value
        if end_index_value < scaled_stop_value:
            print "Specified end of comparison is greater than last value in the data set."
            print "Comparison Ends at Index #:", str(end_index), "with a value of: ", str(end_index_value)
            minmax_end_index = end_index
            break
        else:
            #print "Value 4:",value4
            if value4 < scaled_stop_value:
                rowcounter4 = rowcounter4 + 1
                #print "increment rowcounter4:", rowcounter4
            if value4 >= scaled_stop_value:
                if value4 == scaled_stop_value:
                    row_number4 = rowcounter4
                    print "Comparison Ends at Index #:", str(row_number4), "with a value of: ", str(data_dict[col_name][row_number4])
                    minmax_end_index = row_number4
                    break
                else:
                    row_number4 = rowcounter4 - 1
                    print "Comparison Ends at Index #:", str(row_number4), "with a value of: ", str(data_dict[col_name][row_number4])
                    minmax_end_index = row_number4
                    break
    print "\n*** Start/End Indexes Found ***"
    print "XCol Start Index:", Xcol_start_index
    print "XCol End Index:", Xcol_end_index
    print "minmax Start Index:", minmax_start_index
    print "minmax End Index:", minmax_end_index, "\n"
    return (Xcol_start_index, Xcol_end_index, minmax_start_index, minmax_end_index)

def extract_comp_data(comp_file_info):
    ## Read in d line dict from config file and Process data from source .csv files.
    
    exp_data = []
    mod_data = []
    exp_data_dict = {}
    mod_data_dict = {}
    exp_scatter_data_labels = []
    mod_scatter_data_labels = []
    
    #List of variables from configuration file column names.
    
    exp_data_filename = comp_file_info['Exp_Filename'] #String of filename
    exp_column_name_row_index = int(comp_file_info['Exp_Col_Name_Row'])-1 #Experimental Data Column Name Row Number
    exp_data_row_index = int(comp_file_info['Exp_Data_Row'])-1 #Experimental Data Starting Row Number
    exp_start_data_val = comp_file_info['Exp_Start_(min.)'] #String in minutes to start exp plot data
    exp_stop_data_val = comp_file_info['Exp_End_(min.)'] #String in minutes to stop exp plot data
    exp_start_comp_val = comp_file_info['Exp_Comp_Start_(min.)'] #String in minutes to start exp compare data
    exp_stop_comp_val = comp_file_info['Exp_Comp_End_(min.)'] #String in minutes to start exp compare data
    exp_initial_value = comp_file_info['Exp_Intitial_Value'] #Initial Value for Quantity
    exp_X_column_name_value = comp_file_info['Exp_X_Col_Name'].strip() #Experimental Data X Column Name
    exp_Y_column_name_value = comp_file_info['Exp_Y_Col_Name'].strip() #Experimental Data Y Column Name
    X_Scale_Factor = int(comp_file_info['Scale_X'])
    Y_Scale_Factor = int(comp_file_info['Scale_Y'])
        
    mod_data_filename = comp_file_info['Mod_Filename'] #String of filename
    mod_column_name_row_index = int(comp_file_info['Mod_Col_Name_Row'])-1 #Modeling Data Column Name Row Number
    mod_data_row_index = int(comp_file_info['Mod_Data_Row'])-1 #Modeling Data Starting Row Number
    mod_start_data_val = comp_file_info['Mod_Start_(min.)'] #String in minutes to start mod plot data
    mod_stop_data_val = comp_file_info['Mod_End_(min.)']  #String in minutes to stop mod plot data
    mod_start_comp_val = comp_file_info['Mod_Comp_Start_(min.)'] #String in minutes to start mod compare data
    mod_stop_comp_val = comp_file_info['Mod_Comp_End_(min.)']  #String in minutes to start mod compare data
    mod_initial_value = comp_file_info['Mod_Intitial_Value']       #Initial Value for Quantity
    mod_X_column_name_value = comp_file_info['Mod_X_Col_Name'].strip() #Modeling Data X Column Name
    mod_Y_column_name_value = comp_file_info['Mod_Y_Col_Name'].strip() #Modeling Data Y Column Name
    
    # Create Scatter Data Labels for the comparison results.
    
    if exp_Y_column_name_value[0] == '[':
        print "Exp Column Name List Detected"
        exp_compound_col_names = eval(exp_Y_column_name_value)
        #print "Exp Compound Column Names:", exp_compound_col_names
        for name in exp_compound_col_names:
            print "Exp Sub-Column Name:", name
            exp_scatter_data_labels.append(comp_file_info['Quantity']+"~"+comp_file_info['Group']+"~"+comp_file_info['Dataname']+"~"+name)
    else:
        print "Single Exp. Column Name:", exp_Y_column_name_value
        exp_scatter_data_labels.append(comp_file_info['Quantity']+"~"+comp_file_info['Group']+"~"+comp_file_info['Dataname']+"~"+exp_Y_column_name_value)        
    
    if mod_Y_column_name_value[0] == '[':
        print "Mod Column Name List Detected"
        mod_compound_col_names = eval(mod_Y_column_name_value)
        #print "Mod Compound Column Names:", mod_Y_column_name_value
        for name in mod_compound_col_names:
            print "Mod Sub-Column Name:", name
            mod_scatter_data_labels.append(comp_file_info['Quantity']+"~"+comp_file_info['Group']+"~"+comp_file_info['Dataname']+"~"+name)
    else:
        print "Single Mod. Column Name:", mod_Y_column_name_value
        mod_scatter_data_labels.append(comp_file_info['Quantity']+"~"+comp_file_info['Group']+"~"+comp_file_info['Dataname']+"~"+mod_Y_column_name_value)
    
    #print "Exp Data Labels:\n", exp_scatter_data_labels
    #print "Mod Data Labels:\n", mod_scatter_data_labels
    
    combined_scatter_data_labels = [exp_scatter_data_labels,mod_scatter_data_labels]
    #print "Combined Scatter Data:",combined_scatter_data_labels
    
    min_max = comp_file_info['max/min'] #String indicating if min or max value is required.
    
    group_value = int(comp_file_info['Group'])
    
    try:
        exp_file_object = open(data_directory+exp_data_filename, "U")
    except:
        print "!!! Experimental "+exp_data_filename+" Data File will not open. !!!"
        exit()
        
    try:
        mod_file_object = open(data_directory+mod_data_filename, "U")
    except:
        print "!!! Modeling "+mod_data_filename+" Data File will not open. !!!"
        exit()
    
    ## Start File Processing
    
    #Read in experimental data and flip lists from rows to columns.
    print "Reading in:", exp_data_filename
    exp_data_cols = zip(*csv.reader(exp_file_object))
    #print "exp_data_cols: ",exp_data_cols
    
    #Find X_Axis index number and confirm Col_Name based on Exp_X_Col_Name value in config file.
    column_counter = 0
    for column in exp_data_cols:
        if column[exp_column_name_row_index].strip() == exp_X_column_name_value:
            print "Exp. X Col name is: ",column[exp_column_name_row_index].strip()
            #print "The Index Value is:",column_counter
            exp_Xaxis_column_name = column[exp_column_name_row_index].strip()
        else:
            column_counter = column_counter + 1
            if column_counter == len(exp_data_cols):
                print "!!! Problem with Exp_X_Col_Name value in Config File !!!"
                exit()
            
    #Convert tuples to lists.
    exp_data_list = [list(sublist) for sublist in exp_data_cols]
    
    #Build Experimental Data Dictionary. 
    #Catch errors if conversion of data from string to float fails.
    for exp_list in exp_data_list:
        try:
            temp_list = []
            for x in exp_list[exp_data_row_index:]:
                if x == 'Null' or x == '' or x == 'NaN' or x == 'inf' or x == '-inf':
                    list_value = 'Null'
                else:
                    list_value = float(x)
                temp_list.append(list_value)
            exp_data_dict[exp_list[exp_column_name_row_index].strip()] = temp_list
        except:
            print "!!! Exp Data Conversion in Column Name "+exp_list[exp_column_name_row_index].strip()+". !!!"
            exit()
    
    #Read in model data and flip lists from rows to columns.
    print "Reading in:", mod_data_filename
    mod_data_cols = zip(*csv.reader(mod_file_object))
    #print "mod_data_cols: ",mod_data_cols
    
    #Find X_Axis index number and confirm Col_Name based on Mod_X_Col_Name value in config file.
    column_counter = 0
    for column in mod_data_cols:
        if column[mod_column_name_row_index].strip() == mod_X_column_name_value:
            print "Mod. X Col name is: ",column[mod_column_name_row_index].strip()
            #print "The Index Value is:",column_counter
            mod_Xaxis_column_name = column[mod_column_name_row_index].strip()
        else:
            column_counter = column_counter + 1
            if column_counter == len(mod_data_cols):
                print "!!! Problem with Mod_X_Col_Name value in Config File !!!"
                exit()
    
    #Convert tuples to lists.
    mod_data_list = [list(sublist) for sublist in mod_data_cols]
    
    #Build Prediction/Model Data Dictionary
    #Catch errors if conversion of data from string to float fails.
    for mod_list in mod_data_list:
        try:
            temp_list = []
            for x in mod_list[mod_data_row_index:]:
                if x == 'Null' or x == '' or x == 'NaN' or x == 'inf' or x == '-inf':
                    list_value = 'Null'
                else:
                    list_value = float(x)
                temp_list.append(list_value)
            mod_data_dict[mod_list[mod_column_name_row_index].strip()] = temp_list
        except:
            print "!!! Mod Data Conversion in Column Name "+mod_list[mod_column_name_row_index].strip()+". !!!"
            exit()
    
    # Passing in the X_Axis Column Name.
    exp_comp_ranges = find_start_stop_index(exp_data_dict,exp_Xaxis_column_name,exp_start_data_val,exp_stop_data_val,exp_start_comp_val,exp_stop_comp_val,X_Scale_Factor)
    mod_comp_ranges = find_start_stop_index(mod_data_dict,mod_Xaxis_column_name,mod_start_data_val,mod_stop_data_val,mod_start_comp_val,mod_stop_comp_val,X_Scale_Factor)
    #print "EXP COMP RANGES",exp_comp_ranges
    #print "MOD COMP RANGES",mod_comp_ranges
    
    #### Begin Column specific operations.
    scatter_counter = 0
    
    for scatter_label in combined_scatter_data_labels[0]:
        
        #print scatter_counter
        
        exp_label_temp = []
        mod_label_temp = []
        
        exp_label_temp = split("~",combined_scatter_data_labels[0][scatter_counter])
        mod_label_temp = split("~",combined_scatter_data_labels[1][scatter_counter])
        
        #print "Exp. Label Split:", exp_label_temp
        #print "Mod. Label Split:", mod_label_temp
        
        ##Find max or min values.
        exp_data_values_comp = exp_data_dict[exp_label_temp[3]][exp_comp_ranges[2]:(exp_comp_ranges[3]+1)]
        mod_data_values_comp = mod_data_dict[mod_label_temp[3]][mod_comp_ranges[2]:(mod_comp_ranges[3]+1)]
        
        #print "Exp data values:", exp_data_values_comp
        #print "Mod data values:", mod_data_values_comp
        
        # This allows the d line Quantity value to be set to 0 when either model or experimental data is missing.
        if comp_file_info['Quantity'] == str(0):
            print "Quantity set to 0, no comparison made."
        else:
            if  min_max == 'max':
                print "*** Compute Rise ***"
                temp_exp_data_values = [x for x in exp_data_values_comp if x != 'Null']
                exp_rise_value = max(temp_exp_data_values) - float(exp_initial_value)
                temp_mod_data_values = [x for x in mod_data_values_comp if x != 'Null']
                mod_rise_value = max(temp_mod_data_values) - float(mod_initial_value)
                print "Experimental Initial Value is:", exp_initial_value
                print "Experimental Rise Value is:", exp_rise_value
                print "Model Initial Value is:", mod_initial_value
                print "Model Rise Value is:", mod_rise_value
                print "\n*** Computing Rise Relative Difference ***"
                try:
                    relative_difference = ((mod_rise_value-exp_rise_value)/exp_rise_value)
                    print "Rise Relative Difference is:", relative_difference
                except:
                    print "!!! Computation of Rise relative_difference failed. !!!\nCheck source data for columns listed above."
                    exit()
                    #Append Rise Values to Global Scatter Data Dictionary.
                #print combined_scatter_data_labels[0][scatter_counter]
                scatter_data_dict[combined_scatter_data_labels[0][scatter_counter]] = [exp_rise_value,mod_rise_value,relative_difference]
            elif min_max == 'min':
                print "*** Compute Drop ***"
                temp_exp_data_values = [x for x in exp_data_values_comp if x != 'Null']
                exp_drop_value = float(exp_initial_value) - min(temp_exp_data_values)
                temp_mod_data_values = [x for x in mod_data_values_comp if x != 'Null']
                mod_drop_value = float(mod_initial_value) - min(temp_mod_data_values)
                print "Experimental Initial Value is:", exp_initial_value
                print "Experimental Drop Value is:", exp_drop_value
                print "Model Initial Value is:", mod_initial_value
                print "Model Drop Value is:", mod_drop_value
                print "\n*** Computing Drop Relative Difference ***"
                try:
                    relative_difference = ((mod_drop_value-exp_drop_value)/exp_drop_value)
                    print "Min Relative Difference is:", relative_difference
                except:
                    print "!!! Computation of Min relative_difference failed. !!!\nCheck source data for columns listed above."
                    exit()
                #Append Drop Values to Global Scatter Data Dictionary.
                scatter_data_dict[combined_scatter_data_labels[0][scatter_counter]] = [exp_drop_value,mod_drop_value,relative_difference]
            else:
                print "!!! Min or Max is undefined in the input file. !!!"
                exit()
                
        #Create data lists based on specified ranges
        exp_data_seconds = zip(exp_data_dict[exp_Xaxis_column_name][exp_comp_ranges[0]:(exp_comp_ranges[1]+1)], exp_data_dict[exp_label_temp[3]][exp_comp_ranges[0]:(exp_comp_ranges[1]+1)])
        #print "Exp_data_seconds", exp_data_seconds
        mod_data_seconds = zip(mod_data_dict[mod_Xaxis_column_name][mod_comp_ranges[0]:(mod_comp_ranges[1]+1)], mod_data_dict[mod_label_temp[3]][mod_comp_ranges[0]:(mod_comp_ranges[1]+1)])
        #print "Mod_data_seconds", mod_data_seconds
        
        #Convert time to minutes from seconds.
        exp_data.append([[x[0] / X_Scale_Factor, x[1]] for x in exp_data_seconds])
        #print exp_data
        mod_data.append([[x[0] / X_Scale_Factor, x[1]] for x in mod_data_seconds])
        #print mod_data
        
        #print "\nUpdating the Scatter Counter.\n"
        scatter_counter = scatter_counter + 1
        
    # Close files
    exp_file_object.close()
    mod_file_object.close()
    
    return [exp_data,mod_data]

def comparison_plot(plot_data,exp_data,mod_data):
    #plot_data is a list of values from the 'd' row of the config file being processed.
    
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
    
    # Create line styles that have predetermined color order for each pair in series.  
    # All Experimental data is plotted with solid lines while Model data is dotted.
    expPlotStyle = graph.style.line(lineattrs=[attr.changelist([color.cmyk.Black, color.cmyk.Red, color.cmyk.Green, color.cmyk.Blue]), style.linestyle.solid, style.linewidth(0.06*unit.w_cm)])
    modPlotStyle = graph.style.line(lineattrs=[attr.changelist([color.cmyk.Grey, color.cmyk.Red, color.cmyk.Green, color.cmyk.Blue]), style.linestyle.solid, style.linewidth(0.03*unit.w_cm)])
    
    #Loop strcuture to process compound colum names in d line.
    #print "Exp. Data:",exp_data
    if len(exp_data) > 1 :
        #Set plot legend key text.
        exp_key_list = eval(plot_data['Exp_Key'])
        mod_key_list = eval(plot_data['Mod_Key'])
        exp_plot_counter = 0
        mod_plot_counter = 0
        
        # Loop through and plot Experimental data
        for exp_data_item in exp_data:
            g.plot(graph.data.points(exp_data_item, title=exp_key_list[exp_plot_counter], x=1, y=2),
                  [expPlotStyle])
            exp_plot_counter = exp_plot_counter + 1
            
        # Loop through and plot Experimental data
        for mod_data_item in mod_data:
            g.plot(graph.data.points(mod_data_item, title=mod_key_list[mod_plot_counter], x=1, y=2),
                  [modPlotStyle])
            mod_plot_counter = mod_plot_counter + 1
    else:
        #Set plot legend key text.
        exp_key = plot_data['Exp_Key']
        mod_key = plot_data['Mod_Key']
        
        #print "Exp. Data to Plot:", exp_data[0]
        #print "Mod. Data to Plot:", mod_data[0]
        
        # Plot Experimental data
        g.plot(graph.data.points(exp_data[0], title=exp_key, x=1, y=2),
            [expPlotStyle])
        # Plot Predicted/Model data
        g.plot(graph.data.points(mod_data[0], title=mod_key, x=1, y=2),
            [modPlotStyle])
    
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
    print "\n*** Comparison Plot to: ***\n", plot_file_path+".PDF"

def scatter_plot(group_info,scatter_info,data_set):
    #data_set is a dictionary keyed by quantity, containing lists of groups and X and Y data points.
    
    for quantity_number in scatter_info:
        #print "Dataset for quantity number "+str(quantity_number)+": ", data_set[quantity_number]
        
        if data_set[quantity_number] == []:
            print "No Scatter Plot Data in Quantity "+str(quantity_number)+" Dataset.\n"
        else:
            print "Scatter Plot Data for Quantity "+str(quantity_number)+" Dataset."
            # Set variables for Plot extracted from the first group of lines in config file starting with 'q'.
            
            # Variables for plot.
            plot_title = scatter_info[int(quantity_number)]['Scatter_Plot_Title']
            print plot_title
            
            x_title = scatter_info[int(quantity_number)]['X_Title']
            y_title = scatter_info[int(quantity_number)]['Y_Title']
            min_x = float(scatter_info[int(quantity_number)]['Plot_Min'])
            #print min_x
            max_x = float(scatter_info[int(quantity_number)]['Plot_Max'])
            #print max_x
            min_y = float(scatter_info[int(quantity_number)]['Plot_Min'])
            max_y = float(scatter_info[int(quantity_number)]['Plot_Max'])
            percent_error = float(scatter_info[int(quantity_number)]['%error'])
            title_quadrant = int(scatter_info[int(quantity_number)]['Title_Quadrant'])
            plot_width = int(scatter_info[int(quantity_number)]['Plot_Width(cm)'])
            # Specify the position and line spacing of the plot key.
            key_pos = scatter_info[int(quantity_number)]['Key_Position']
            key_dist = 0.15*unit.v_cm
            
            #Create filename from fields in input file record.
            plot_file_name = scatter_info[int(quantity_number)]['Plot_Filename']
            #print plot_file_name
            
            # Determine the location for the key, alignment based on key_quadrant setting.
            if key_pos == "tl" or "tc" or "tr" or "ml" or "mc" or "mr" or "bl" or "bc" or "br":
                ()
                #print "Key Position =", key_pos
            else:
                print "The key position was not specified.\nUsing the default bottom right position."
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
            
            if min_x < 0:
                #print errorLineCenterPoints
                lower_bound = ((min_y)+((min_y)*(percent_error / 100)))
                #print "Lower Bound:", lower_bound
                errorLineLowerPoints = [[min_x,lower_bound],[max_x,max_y]]
                #print "Lower Error Line Points:", errorLineLowerPoints
                upper_bound = ((min_y)-((min_y)*(percent_error/100)))
                #print "Upper Bound:", 
                errorLineUpperPoints = [[min_x,upper_bound],[max_x,max_y]]
                #print "Upper Error Line Points:", errorLineUpperPoints
            else:
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
                    
            if percent_error == 0:
                print "No Error Bars Drawn"
            else:
                g.plot(graph.data.points(errorLineLowerPoints, title=None, x=1, y=2),
                        [graph.style.line([style.linewidth.Thin, style.linestyle.dashed])])
                    
                g.plot(graph.data.points(errorLineUpperPoints, title=None, x=1, y=2),
                        [graph.style.line([style.linewidth.Thin, style.linestyle.dashed])])
            
            #One point at a time added to plot from each data set.
            # Iterate over items in scatter data dictionary key for items that are not [].
            # Append data sets to scatter_plot_data_list
            # colors for symbols are from http://pyx.sourceforge.net/manual/colorname.html
            
            scatter_plot_data_list = []
            
            grouped_data = {}
            grouped_data_list = range(len(group_quantity_data_dicts[0])+1)
            #print "Grouped Data List:", grouped_data_list
            
            #print "DataSet for Quantity "+str(quantity_number)+":", data_set[quantity_number]
            if len(data_set[quantity_number]) > 1:
                #print "Grouped Scatter Data:"
                #print data_set[quantity_number]
                for arr_temp in grouped_data_list:
                    grouped_data_list[arr_temp] = []
                for data_set_item in data_set[quantity_number]:
                    #print data_set_item
                    #print "Data for group "+data_set_item[0]+":", data_set_item[1]
                    grouped_data_list[int(data_set_item[0])].append(data_set_item[1])
                #print "Grouped data list:", grouped_data_list
                                
                group_counter = 0
                for j in grouped_data_list:
                    #print "J =", j
                    if j != []:
                        #print group_counter
                        
                        # Pull group symbol specifications from config file.
                        config_group_symbol = group_info[group_counter]["Symbol"]
                        #print config_group_symbol
                        group_symbol = "graph.style.symbol."+config_group_symbol
                        #print group_symbol
                        
                        config_group_symbol_color = group_info[group_counter]["Color"]
                        #print config_group_symbol_color
                        #group_color = "color.rgb."+config_group_symbol_color
                        
                        config_group_symbol_filled = group_info[group_counter]["Filled"]
                        #print config_group_symbol_filled
                        
                        if config_group_symbol_filled == 'yes':
                            fillstyle = "deco.filled([color.cmyk."+config_group_symbol_color+"])"
                        else:
                            fillstyle = "deco.stroked([color.cmyk."+config_group_symbol_color+"])"
                        #print group_symbol, fillstyle     
                        
                        #Create temporary symbol style.
                        tempstyle = "graph.style.symbol("+group_symbol+", size=0.1*unit.v_cm, symbolattrs=["+fillstyle+"])"
                        #print "TempStyle:", tempstyle
                        
                        scatterpointstyle = eval(tempstyle)
                        
                        #print group_info[group_counter]["Group_Title"]
                        g.plot(graph.data.points(j, x=1, y=2, title=group_info[group_counter]["Group_Title"]), [scatterpointstyle])
                    else:
                        pass
                    group_counter = group_counter + 1
                    
            else:
                print "Non-Grouped Scatter Data:"
                #print data_set[quantity_number]
                scatter_plot_data = []
                scatter_plot_data.append(data_set[quantity_number][0][1])
                #print scatter_plot_data
                
            #print grouped_data_list
            
            # Now plot the Title text, alignment based on title quadrant setting.
            if title_quadrant == 1:
                g.text(0.1, g.height - 0.2, plot_title, [text.halign.left, text.valign.top, text.size.small])
            elif title_quadrant == 2:
                g.text(g.width-0.1, g.height - 0.2, plot_title, [text.halign.right, text.valign.top, text.size.small])
            elif title_quadrant == 3:
                g.text(0.1, 0.2, plot_title, [text.halign.left, text.valign.bottom, text.size.small])
            elif title_quadrant == 4:
                g.text(g.width-0.1, 0.2, plot_title, [text.halign.right, text.valign.bottom, text.size.small])
            else:
                print "A title location was not specified./nUsing the default top left quadrant."
                g.text(0.1, g.height - 0.2, plot_title, [text.halign.left, text.valign.top, text.size.small])
            
            #Make %error text on plot by error bars.
            # pos_percent_error = str(percent_error)+"%"
            # neg_percent_error = "-"+str(percent_error)+"%"
            # g.text(g.width - 0.4, g.height - 0.3, pos_percent_error, [text.halign.center, text.valign.middle, text.size.tiny])
            # g.text(g.width - 0.2, g.height - 0.4, neg_percent_error, [text.halign.center, text.valign.middle, text.size.tiny])
            
            # Write the output
            plot_file_path = output_directory+plot_file_name
            #print plot_file_path
            g.writePDFfile(plot_file_path)
            print "Scatter Plot to: \n", plot_file_path+".PDF\n"


### Start of Main Code
print "**** READING CONFIGURATION FILE ****"

##Get information from config file.
group_quantity_data_dicts = extract_config_data(config_file_name)
print "\nThere are "+str(len(group_quantity_data_dicts[0]))+" scatter data groups, (g lines)."
#print group_quantity_data_dicts[0]
print "There are "+str(len(group_quantity_data_dicts[1]))+" scatter data quantities to plot, (q lines)."
#print group_quantity_data_dicts[1]
print "There are "+str(len(group_quantity_data_dicts[2]))+" comparison data sets to plot, (d lines)."

## Create comparison plots
print "**** CREATING COMPARISON PLOTS ****"
d_count = 1
for data_record in group_quantity_data_dicts[2]:
    # Each d line, data_record, may contain compound column names from the config file.
    print "*** #"+str(d_count)+" of "+str(len(group_quantity_data_dicts[2]))+" comparison records. ***\n"
    # Extract relevant portions of comparison data as defined in config file.
    comp_data_to_plot = extract_comp_data(group_quantity_data_dicts[2][data_record])
    #print "Comparison Data to Plot:", comp_data_to_plot
    
    #Seperate experimental and model data lists.
    exp_plot_data = comp_data_to_plot[0]
    mod_plot_data = comp_data_to_plot[1]
    
    #print "Exp Plot Data:", exp_plot_data
    #print "Mod Plot Data:", mod_plot_data
    
    # Create plot for data_record.
    comparison_plot(group_quantity_data_dicts[2][data_record],exp_plot_data,mod_plot_data)
    d_count = d_count + 1
    print "\n"

## Create scatter plots
print "**** CREATING SCATTER PLOTS ****"
scatter_quantity = 1
scatter_group = 1
temp_scatter_data_list = []

#Grouping Scatter Data by Quantity
for scatter_plot_record in sorted(group_quantity_data_dicts[1]):
    #print "Scatter Plot:", scatter_plot_record
    #print "Quantity:", group_quantity_data_dicts[1][scatter_plot_record]['Comparison_Quantities']
        
    for scatter_data_key in sorted(scatter_data_dict):
        split_key = split("~",scatter_data_key)
        if split_key[0] == str(scatter_plot_record) and scatter_data_dict[scatter_data_key] != []:
            temp_scatter_data_list.append([split_key[1], scatter_data_dict[scatter_data_key][:2]])
        else:
            pass
    
    #print temp_scatter_data_list
    combined_scatter_data[scatter_quantity] = temp_scatter_data_list
    temp_scatter_data_list = []
    scatter_quantity = scatter_quantity + 1

# Plot Data
#print "Data to Scatter Plot:", combined_scatter_data
scatter_plot(group_quantity_data_dicts[0],group_quantity_data_dicts[1],combined_scatter_data)

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

print "***  Processing finished, thank you for your patience.  ***"
