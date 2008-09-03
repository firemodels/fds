import csv
from re import *
from pyx import *

### Set Global Variables for Data to Process
## Validation Data Set (1), Verification Data Set (2), Examples Data Set (3), Trainier Data Set (4)
data_set = 1 # Input Value

## Validation Data
if data_set == 1:
    data_directory = "../../Validation/"
    output_directory = "../../Manuals/FDS_5_Validation_Guide/FIGURES/"
    config_file_name = "Validation_Data_Config_File.csv"
    print "**** Processing Validation Data Set ****\n"

## Verification Data
if data_set == 2:
    data_directory = "../../Verification/"
    output_directory = "../../Manuals/FDS_5_Verification_Guide/FIGURES/"
    config_file_name = "Verification_Data_Config_File.csv"
    print "**** Processing Verification Data Set ****\n"

## Examples Data (Currently not working because the file system is changed)
if data_set == 3:
    data_directory = "../../Verification/"
    output_directory = "../../Manuals/FDS_5_User_Guide/FIGURES/"
    config_file_name = "Examples_Data_Config_File.csv"
    print "**** Processing Examples Data Set ****\n"

## Trainier Data
if data_set == 4:
    data_directory = "../../Training/"
    output_directory = "../../Manuals/FDS_SMV_Training_Guide/datafigures/"
    config_file_name = "Training_Examples_Data_Config_File.csv"
    print "**** Processing Training Data Set ****\n"


### Set Diagnostic Output Level
## Uncomment the level of diagnostics desired.
## 1 = Minimal, 2 = Normal, 3 = Maximum.
diagnostic_level = 2 # Input Value

print "**** Diagnostics Set at Level", diagnostic_level, "****\n"


### Set What Plots to Create: 
## BOTH (1), Comparison Only (2), Scatter Only (3)
process_set = 1 # Input Value

if diagnostic_level >= 1:
    if process_set == 1:
        print "**** Plotting both Comparison and Scatter Data ****\n"
    if process_set == 2:
        print "**** Plotting Only Comparison Data ****\n"
    if process_set == 3:
        print "**** Plotting Only Scatter Data ****\n"

### Set character to be used for indicating the comparison data set in the config file.
## The default character is 'd' if you change the value below to something else 
## then only the lines starting with that character will be read in for processing.
## NOTE: You must also change the line containing the d line column names to the same character string.
data_line_char = 'd' # Input Value
print "**** Data Character Set to '"+data_line_char+"' ****\n"


### Set Global Variables to Zero or Null Sets.
scatter_data_dict = {}
combined_scatter_data = {}

### Define Functions

def extract_config_data(config_file):
    if diagnostic_level >= 1:
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
    if diagnostic_level >= 2:
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
            if diagnostic_level >= 3:
                print "   <3> This is group line #"+str(group_counter)+"."
        elif list_item[0] == 'q':
            if quantity_counter < 1:
                quantity_header = list_item[2:]
            if quantity_counter >= 1:
                for x in range(len(quantity_header)):
                    keyed_quantities[quantity_header[x]] = list_item[x+2]
                if diagnostic_level >= 3:
                    print "   <3> List item 1:",int(list_item[1])
                quantities_dict[int(list_item[1])] = keyed_quantities
                keyed_quantities = {}
            quantity_counter += 1
        elif list_item[0] == data_line_char:
            if data_counter < quantity_counter:
                data_counter = quantity_counter - 1
                data_header = list_item[1:]
            if data_counter >= quantity_counter:
                for x in range(len(data_header)):
                    keyed_data[data_header[x].strip()] = list_item[x+1]
                data_key_name = keyed_data['Quantity'].strip()+"~"+keyed_data['Group'].strip()+"~"+keyed_data['Dataname'].strip()+"~"+keyed_data['d1_Dep_Col_Name'].strip()+"~"+keyed_data['d1_Ind_Col_Name'].strip()
                if diagnostic_level >= 3:
                    print "   <3> Key Name:", data_key_name
                data_dict[data_key_name] = keyed_data
                if diagnostic_level >= 3:
                    print "   <3> Data for Key Name Above:", data_dict[data_key_name]
                keyed_data = {}
            data_counter += 1
        else:
            skip_counter = skip_counter + 1
            if diagnostic_level >= 3:
                print "   <3> No g, q or "+data_line_char+", skipping row."            
    if diagnostic_level >= 2:
        print "There were "+str(skip_counter)+" lines skipped, out of the "+str(len(config_lists))+" lines read in."
    #Return a single list object containing the dictionaries.
    if diagnostic_level >= 3:
        print "   <3> Groups Dictionary: ", groups_dict
        print "   <3> Quantities Dictionary: ", quantities_dict
        print "   <3> Data Dictionary: ", data_dict
    return [groups_dict,quantities_dict,data_dict]

def find_start_stop_index(data_dict,col_name,start_data,stop_data,start_comp,stop_comp,ind_scale):
    #This function is used to find index numbers for start and stop points in plotting and metric values.
    if diagnostic_level >= 2:
        print "\n*** Find Start/End Indexes for Independent Data Column ***"
    if diagnostic_level >= 3:
        print "   <3> Start Data:", start_data
        print "   <3> Stop Data:", stop_data
        print "   <3> Independent Data Scale Factor:", ind_scale
        
    if diagnostic_level >= 2:
        print "Independent Data Column name:", col_name
    if diagnostic_level >= 3:
        print "   <3> Independent Data Column Data:", data_dict[col_name]
    ## Find Start Value Index Number    
    rowcounter1 = 0
    if diagnostic_level >= 2:
        print "*** Finding Independent Data Column Start Index ***"
    for value1 in data_dict[col_name]:
        if diagnostic_level >= 3:
            print "   <3> Value 1:", value1, " : ", rowcounter1
        if value1 == 'Null':
            continue
        else:
            if value1 >= (float(start_data)*float(ind_scale)):
                if diagnostic_level >= 3:
                    print "   <3> Independent Data Column Starts at row #:", str(rowcounter1)
                    print "   <3> With a value of:", str(value1)
                ind_col_start_index = rowcounter1
                break
        rowcounter1 = rowcounter1 + 1
    ## Find End of Independent Data Column, then Index Value
    rowcounter2 = 0
    end_index = 0
    index_count = 0
    if diagnostic_level >= 2:
        print "*** Finding Index of last good value in Independent Data Column ***"    
    for end_data_index in data_dict[col_name]:
        if diagnostic_level >= 3:
            print "   <3> ", end_data_index, " : ", index_count
        if end_data_index == 'Null':
            end_index = index_count - 1
            break
        else:
            index_count = index_count + 1
    if end_index == 0:
        end_index = index_count - 1
    if diagnostic_level >= 3:
        print "   <3> Length of Independent Data Column:", end_index + 1
        print "   <3> Last value in Independent Data Column:", data_dict[col_name][end_index]
        print "   <3> Index of Last Numeric Value:", end_index
    if diagnostic_level >= 2:
        print "*** Finding Independent Data Column End Index ***"
    for value2 in data_dict[col_name]:
        if diagnostic_level >= 3:
            print "Value 2:", value2, " : ", rowcounter2
        if float(data_dict[col_name][end_index]) < (float(stop_data)*float(ind_scale)):
            if diagnostic_level >= 2:
                print "Specified end of Independent Data Column Data is greater than end of time in the data set. \nUsing last value in the Independent Data Column.\n"
                print "Value used is: "+str(float(data_dict[col_name][end_index]))+"\n"
            ind_col_end_index = (end_index)
            break
        if value2 < (float(stop_data)*float(ind_scale)):
            rowcounter2 = rowcounter2 + 1
        else:
            row_number2 = (rowcounter2) - 1
            if diagnostic_level >= 2:
                print "Independent Data Column Ends at Index #: ", str(row_number2), "with a value of: ", str(data_dict[col_name][row_number2])
            ind_col_end_index = row_number2
            break
    #Find Comparison Start Index for Independent Data
    rowcounter3 = 0
    if diagnostic_level >= 2:
        print "*** Finding Start Index for Metric Data ***"
    for value3 in data_dict[col_name]:
        if diagnostic_level >= 3:
            print "Value 3:", value3, " : ", rowcounter3
        if value3 >= (float(start_comp)*float(ind_scale)):
            if diagnostic_level >= 2:
                print "Metric Data starts at Index #:", str(rowcounter3), "with a value of:", str(value3)
            metric_start_index = rowcounter3
            break
        rowcounter3 = rowcounter3 + 1 
    #Find Comparison End Index for Independent Data
    rowcounter4 = 0
    if diagnostic_level >= 2:
        print "*** Finding End Index for Metric Data ***"
    for value4 in data_dict[col_name]:
        scaled_stop_value = (float(stop_comp)*float(ind_scale))
        end_index_value = float(data_dict[col_name][end_index])
        if diagnostic_level >= 3:
            print "   <3> Scaled Stop Value and End Index Value:", scaled_stop_value, " : ",end_index_value
        if end_index_value < scaled_stop_value:
            if diagnostic_level >= 2:
                print "Specified end of Metric data is greater than last value in the Independent Data Column.\nUsing last value in the Independent Data Column.\n"
                print "Metric Data ends at Index #:", str(end_index), "with a value of: ", str(end_index_value)
            metric_end_index = end_index
            break
        else:
            if diagnostic_level >= 3:
                print "Value 4:", value4, " : ", rowcounter4
            if value4 < scaled_stop_value:
                rowcounter4 = rowcounter4 + 1
            if value4 >= scaled_stop_value:
                if value4 == scaled_stop_value:
                    row_number4 = rowcounter4
                    if diagnostic_level >= 2:
                        print "Comparison Ends at Index #:", str(row_number4), "with a value of: ", str(data_dict[col_name][row_number4])
                    metric_end_index = row_number4
                    break
                else:
                    row_number4 = rowcounter4 - 1
                    if diagnostic_level >= 2:
                        print "Comparison Ends at Index #:", str(row_number4), "with a value of: ", str(data_dict[col_name][row_number4])
                    metric_end_index = row_number4
                    break
    if diagnostic_level >= 2:
        print "\n*** Start/End Indexes Found ***"
        print "Independent Data Col Start Index:", ind_col_start_index
        print "Independent Data Col End Index:", ind_col_end_index
        print "Metric Start Index:", metric_start_index
        print "Metric End Index:", metric_end_index, "\n"
    return (ind_col_start_index, ind_col_end_index, metric_start_index, metric_end_index)

def extract_comp_data(comp_file_info):
    ## Read in d line dict from config file and Process data from source .csv files.
    
    d1_data = []
    d2_data = []
    d1_data_dict = {}
    d2_data_dict = {}
    d1_scatter_data_labels = []
    d2_scatter_data_labels = []
    
    #List of variables from configuration file column names.
    
    d1_data_filename = comp_file_info['d1_Filename'] #String of filename
    d1_column_name_row_index = int(comp_file_info['d1_Col_Name_Row'])-1 #Data 1, Column Name Row Number
    d1_data_row_index = int(comp_file_info['d1_Data_Row'])-1 #Data 1, Starting Row Number
    d1_start_data_val = comp_file_info['d1_Start'] #String value to start d1 plot data
    d1_stop_data_val = comp_file_info['d1_End'] #String value to stop d1 plot data
    d1_start_comp_val = comp_file_info['d1_Comp_Start'] #String value to start d1 compare data
    d1_stop_comp_val = comp_file_info['d1_Comp_End'] #String value to start d1 compare data
    d1_initial_value = comp_file_info['d1_Initial_Value'] #Initial Value for Quantity
    d1_ind_column_name_value = comp_file_info['d1_Ind_Col_Name'].strip() #Data 1, Independent Data Column Name
    d1_Dep_column_name_value = comp_file_info['d1_Dep_Col_Name'].strip() #Data 1, Dep Column Name
    ind_Scale_Factor = float(comp_file_info['Scale_Ind'])
    Dep_Scale_Factor = float(comp_file_info['Scale_Dep'])
        
    d2_data_filename = comp_file_info['d2_Filename'] #String of filename
    d2_column_name_row_index = int(comp_file_info['d2_Col_Name_Row'])-1 #Data Set 2, Data Column Name Row Number
    d2_data_row_index = int(comp_file_info['d2_Data_Row'])-1 #Data Set 2, Data Starting Row Number
    d2_start_data_val = comp_file_info['d2_Start'] #String value to start d2 plot data
    d2_stop_data_val = comp_file_info['d2_End']  #String value to stop d2 plot data
    d2_start_comp_val = comp_file_info['d2_Comp_Start'] #String value to start d2 compare data
    d2_stop_comp_val = comp_file_info['d2_Comp_End']  #String value to start d2 compare data
    d2_initial_value = comp_file_info['d2_Initial_Value']  #Initial value for Quantity
    d2_ind_column_name_value = comp_file_info['d2_Ind_Col_Name'].strip() #Data Set 2, Independent Data Column Name
    d2_Dep_column_name_value = comp_file_info['d2_Dep_Col_Name'].strip() #Data Set 2, Dep Column Name
    
    # Create Scatter Data Labels for the comparison results.
    
    if d1_Dep_column_name_value[0] == '[':
        if diagnostic_level >= 2:
            print "Data Set 1, Column Name List Detected"
        d1_compound_col_names = eval(d1_Dep_column_name_value)
        if diagnostic_level >= 3:
            print "   <3> Data Set 1, Compound Column Names:", d1_compound_col_names
        for name in d1_compound_col_names:
            if diagnostic_level >= 2:
                print "Data Set 1, Sub-Column Name:", name
            d1_scatter_data_labels.append(comp_file_info['Quantity']+"~"+comp_file_info['Group']+"~"+comp_file_info['Dataname']+"~"+name)
    else:
        if diagnostic_level >= 2:
            print "Single Data Set 1, Column Name:", d1_Dep_column_name_value
        d1_scatter_data_labels.append(comp_file_info['Quantity']+"~"+comp_file_info['Group']+"~"+comp_file_info['Dataname']+"~"+d1_Dep_column_name_value+"~"+d1_ind_column_name_value)        
    
    if d2_Dep_column_name_value[0] == '[':
        if diagnostic_level >= 2:
            print "Data Set 2, Column Name List Detected"
        d2_compound_col_names = eval(d2_Dep_column_name_value)
        if diagnostic_level >= 3:
            print "   <3> Data Set 2, Compound Column Names:", d2_Dep_column_name_value
        for name in d2_compound_col_names:
            if diagnostic_level >= 2:
                print "Data Set 2, Sub-Column Name:", name
            d2_scatter_data_labels.append(comp_file_info['Quantity']+"~"+comp_file_info['Group']+"~"+comp_file_info['Dataname']+"~"+name)
    else:
        if diagnostic_level >= 2:
            print "Single Data Set 2, Column Name:", d2_Dep_column_name_value
        d2_scatter_data_labels.append(comp_file_info['Quantity']+"~"+comp_file_info['Group']+"~"+comp_file_info['Dataname']+"~"+d2_Dep_column_name_value+"~"+d2_ind_column_name_value)
    
    if diagnostic_level >= 3:
        print "   <3> Data Set 1, Data Labels:\n", d1_scatter_data_labels
        print "   <3> Data Set 2, Data Labels:\n", d2_scatter_data_labels
    
    combined_scatter_data_labels = [d1_scatter_data_labels,d2_scatter_data_labels]
    if diagnostic_level >= 3:
        print "   <3> Combined Scatter Data:",combined_scatter_data_labels
    
    metric = comp_file_info['Metric'] #String indicating the type of metric required.
    
    group_value = int(comp_file_info['Group'])
    
    try:
        d1_file_object = open(data_directory+d1_data_filename, "U")
    except:
        print "!!! Data Set 1, filename "+d1_data_filename+" will not open. !!!"
        exit()
        
    try:
        d2_file_object = open(data_directory+d2_data_filename, "U")
    except:
        print "!!! Data Set 2, filename "+d2_data_filename+" will not open. !!!"
        exit()
    
    if diagnostic_level >= 2:
        print "*** Start File Processing ***"
    
    # Read in Data Set 1, data and flip lists from rows to columns.
    if diagnostic_level >= 2:
        print "Reading in:", d1_data_filename
    for x in range(d1_column_name_row_index):
        d1_file_object.next()
    d1_data_cols = zip(*csv.reader(d1_file_object))
    if diagnostic_level >= 3:
        print "   <3> Data Set 1, Data Columns:",d1_data_cols
    
    # Find Ind_Axis index number and confirm Col_Name based on d1_Ind_Col_Name value in config file.
    column_counter = 0
    for column in d1_data_cols:
        if column[0].strip() == d1_ind_column_name_value:
            if diagnostic_level >= 2:
                print "Data Set 1, Independent Data Col name is: ",column[0].strip()
            if diagnostic_level >= 3:
                print "   <3> The Index Value is:",column_counter
            d1_ind_axis_column_name = column[0].strip()
        else:
            column_counter = column_counter + 1
            if column_counter == len(d1_data_cols):
                print "!!! Problem with d1_Ind_Col_Name: "+column[0].strip()+" value in Config File !!!"
                exit()
            
    #Convert tuples to lists.
    d1_data_list = [list(sublist) for sublist in d1_data_cols]
    if diagnostic_level >= 3:
        print "   <3> Data Set 1, Data List:", d1_data_list
    
    if diagnostic_level >= 2:
        print "*** Build Data Set 1 Dictionary. ***" 
    #Catch errors if conversion of data from string to float fails.
    for d1_list in d1_data_list:
        if diagnostic_level >= 3:
            print "   <3> Data Set 1, List:", d1_list
        try:
            temp_list = []
            for x in d1_list[(d1_data_row_index-d1_column_name_row_index):]:
                if x == 'Null' or x == '' or x == 'NaN' or x == 'inf' or x == '-inf':
                    list_value = 'Null'
                else:
                    list_value = float(x)
                temp_list.append(list_value)
                if diagnostic_level >= 3:
                    print "   <3> Temp List:", temp_list
            d1_data_dict[d1_list[0].strip()] = temp_list
        except:
            print "!!! Data Set 1, Conversion Error in Column Name "+d1_list[0].strip()+". !!!"
            exit()
    
    #Read in d2 data and flip lists from rows to columns.
    if diagnostic_level >= 2:
        print "Reading in:", d2_data_filename
    for x in range(d2_column_name_row_index):
        d2_file_object.next()
    d2_data_cols = zip(*csv.reader(d2_file_object))
    if diagnostic_level >= 3:
        print "   <3> Data Set 2, Data Columns:", d2_data_cols
        
    #Find Ind_Axis index number and confirm Col_Name based on d2_Ind_Col_Name value in config file.
    column_counter = 0
    for column in d2_data_cols:
        if column[0].strip() == d2_ind_column_name_value:
            if diagnostic_level >= 2:
                print "Data Set 2, Independent Data Col name is: ",column[0].strip()
            if diagnostic_level >= 3:
                print "   <3> The Index Value is: ",column_counter
            d2_ind_axis_column_name = column[0].strip()
        else:
            column_counter = column_counter + 1
            if column_counter == len(d2_data_cols):
                print "!!! Problem with d2_Ind_Col_Name value in Config File !!!"
                exit()
    
    #Convert tuples to lists.
    d2_data_list = [list(sublist) for sublist in d2_data_cols]
    
    #Build Prediction/Data Set 2, Data Dictionary
    #Catch errors if conversion of data from string to float fails.
    for d2_list in d2_data_list:
        try:
            temp_list = []
            for x in d2_list[(d2_data_row_index-d2_column_name_row_index):]:
                if x == 'Null' or x == '' or x == 'NaN' or x == 'inf' or x == '-inf':
                    list_value = 'Null'
                else:
                    list_value = float(x)
                temp_list.append(list_value)
            d2_data_dict[d2_list[0].strip()] = temp_list
        except:
            print "!!! Data Set 2, Conversion Error in Column Name "+d2_list[0].strip()+". !!!"
            exit()
    
    # Passing in the Ind_Axis Column Name.
    d1_comp_ranges = find_start_stop_index(d1_data_dict,d1_ind_axis_column_name,d1_start_data_val,d1_stop_data_val,d1_start_comp_val,d1_stop_comp_val,ind_Scale_Factor)
    d2_comp_ranges = find_start_stop_index(d2_data_dict,d2_ind_axis_column_name,d2_start_data_val,d2_stop_data_val,d2_start_comp_val,d2_stop_comp_val,ind_Scale_Factor)
    if diagnostic_level >= 3:
        print "   <3> D1 COMP RANGES: ",d1_comp_ranges
        print "   <3> D2 COMP RANGES: ",d2_comp_ranges
    
    #### Begin Column specific operations.
    scatter_counter = 0
    
    for scatter_label in combined_scatter_data_labels[0]:
        
        if diagnostic_level >= 3:
            print "   <3> Scatter Counter Value:", scatter_counter
        
        d1_label_temp = []
        d2_label_temp = []
        
        d1_label_temp = split("~",combined_scatter_data_labels[0][scatter_counter])
        d2_label_temp = split("~",combined_scatter_data_labels[1][scatter_counter])
        
        if diagnostic_level >= 3:
            print "   <3> Data Set 1, Label Split:", d1_label_temp
            print "   <3> Data Set 2, Label Split:", d2_label_temp
        
        ##Find metric values.
        d1_data_values_comp = d1_data_dict[d1_label_temp[3]][d1_comp_ranges[2]:(d1_comp_ranges[3]+1)]
        d2_data_values_comp = d2_data_dict[d2_label_temp[3]][d2_comp_ranges[2]:(d2_comp_ranges[3]+1)]
        
        if diagnostic_level >= 3:
            print "   <3> Data Set 1, data values:", d1_data_values_comp
            print "   <3> Data Set 2, data values:", d2_data_values_comp
        
        # This allows the d line Quantity value to be set to 0 when either d1 or d2 data is missing.
        if comp_file_info['Quantity'] == str(0):
            print "Quantity set to 0, no comparison made."
        else:
            if  metric == 'max':
                if diagnostic_level >= 2:
                    print "*** Compute Rise ***"
                temp_d1_data_values = [x for x in d1_data_values_comp if x != 'Null']
                d1_rise_value = max(temp_d1_data_values) - float(d1_initial_value)
                temp_d2_data_values = [x for x in d2_data_values_comp if x != 'Null']
                d2_rise_value = max(temp_d2_data_values) - float(d2_initial_value)
                if diagnostic_level >= 2:
                    print "Data Set 1, Initial Value is:", d1_initial_value
                    print "Data Set 1, Rise Value is:", d1_rise_value
                    print "Data Set 2, Initial Value is:", d2_initial_value
                    print "Data Set 2, Rise Value is:", d2_rise_value
                    print "\n*** Computing Rise Relative Difference ***"
                try:
                    relative_difference = ((d2_rise_value-d1_rise_value)/d1_rise_value)
                    if diagnostic_level >= 2:
                        print "Rise Relative Difference is:", relative_difference
                except:
                    print "!!! Computation of Rise relative_difference failed. !!!\nCheck source data for columns listed above."
                    exit()
                    #Append Rise Values to Global Scatter Data Dictionary.
                if diagnostic_level >= 3:
                    print "   <3> Scatter Data Labels:", combined_scatter_data_labels[0][scatter_counter]
                scatter_data_dict[combined_scatter_data_labels[0][scatter_counter]] = [d1_rise_value,d2_rise_value,relative_difference]
            elif metric == 'min':
                if diagnostic_level >= 2:
                    print "*** Compute Drop ***"
                temp_d1_data_values = [x for x in d1_data_values_comp if x != 'Null']
                d1_drop_value = float(d1_initial_value) - min(temp_d1_data_values)
                temp_d2_data_values = [x for x in d2_data_values_comp if x != 'Null']
                d2_drop_value = float(d2_initial_value) - min(temp_d2_data_values)
                if diagnostic_level >= 2:
                    print "Data Set 1, Initial Value is:", d1_initial_value
                    print "Data Set 1, Drop Value is:", d1_drop_value
                    print "Data Set 2, Initial Value is:", d2_initial_value
                    print "Data Set 2, Drop Value is:", d2_drop_value
                    print "\n*** Computing Drop Relative Difference ***"
                try:
                    relative_difference = ((d2_drop_value-d1_drop_value)/d1_drop_value)
                    if diagnostic_level >= 2:
                        print "Min Relative Difference is:", relative_difference
                except:
                    print "!!! Computation of Min relative_difference failed. !!!\nCheck source data for columns listed above."
                    exit()
                #Append Drop Values to Global Scatter Data Dictionary.
                scatter_data_dict[combined_scatter_data_labels[0][scatter_counter]] = [d1_drop_value,d2_drop_value,relative_difference]
            else:
                print "!!! Metric is undefined in the input file. !!!"
                exit()
                
        #Create data lists based on specified ranges
        d1_data_seconds = zip(d1_data_dict[d1_ind_axis_column_name][d1_comp_ranges[0]:(d1_comp_ranges[1]+1)], d1_data_dict[d1_label_temp[3]][d1_comp_ranges[0]:(d1_comp_ranges[1]+1)])
        if diagnostic_level >= 3:
            print "   <3> Data Set 1, Data:", d1_data_seconds
        d2_data_seconds = zip(d2_data_dict[d2_ind_axis_column_name][d2_comp_ranges[0]:(d2_comp_ranges[1]+1)], d2_data_dict[d2_label_temp[3]][d2_comp_ranges[0]:(d2_comp_ranges[1]+1)])
        if diagnostic_level >= 3:
            print "   <3> Data Set 2, Data:", d2_data_seconds
        
        #Scale Ind_Axis Data.
        d1_data.append([[x[0] / ind_Scale_Factor, x[1]] for x in d1_data_seconds])
        if diagnostic_level >= 3:
            print "   <3> Scaled Data Set 1, Data:", d1_data
        d2_data.append([[x[0] / ind_Scale_Factor, x[1]] for x in d2_data_seconds])
        if diagnostic_level >= 3:
            print "   <3> Scaled Prediction Data:", d2_data
            
        #Need to Scale Dep_Axis Data...
        
        scatter_counter = scatter_counter + 1
        if diagnostic_level >= 3:
            print "\n   <3> Scatter Counter:", scatter_counter, "\n"
    
    # Close files
    d1_file_object.close()
    d2_file_object.close()
    
    return [d1_data,d2_data]

def comparison_plot(plot_data,d1_data,d2_data):
    #plot_data is a list of values from the 'd' row of the config file being processed.
    
    # Variables for plot.
    plot_title = plot_data['Plot_Title']
    if diagnostic_level >= 3:
        print "   <3> Plot Title:", plot_title
    ind_title = plot_data['Ind_Title']
    dep_title = plot_data['Dep_Title']   
    min_ind = float(plot_data['Min_Ind'])
    max_ind = float(plot_data['Max_Ind'])
    min_dep = float(plot_data['Min_Dep'])
    max_dep = float(plot_data['Max_Dep'])
    title_position = plot_data['Title_Position']
    key_pos = plot_data['Key_Position']
    key_dist = 0.2*unit.v_cm
    v_dist   = 0.7*unit.v_cm
    h_dist   = 0.4*unit.v_cm
    plot_width = float(plot_data['Plot_Width(cm)'])
    
    #Create filename from fields in input file record.
    plot_file_name = plot_data["Plot_Filename"]
    
    # Determine the location for the key, alignment based on key_pos setting.
    if key_pos == "tl" or "tc" or "tr" or "ml" or "mc" or "mr" or "bl" or "bc" or "br":
        ()
        if diagnostic_level >= 3:
            print "   <3> Key Position =", key_pos
    else:
        print "The key position was not specified.\nUsing the default bottom right position."
        key_pos = "br"
    
    # Determine the location for the title, alignment based on Title_Position setting.
    if title_position == "tl" or "tc" or "tr" or "ml" or "mc" or "mr" or "bl" or "bc" or "br":
        ()
        if diagnostic_level >= 3:
            print "   <3> Title Position =", title_position
    else:
        if diagnostic_level >= 2:
            print "The Title position was not specified.\nUsing the default top left position."
        key_pos = "tl"

    #Begin Plotting
    # Initialize graph object
    g = graph.graphxy(width=plot_width, ratio=4./3, key=graph.key.key(pos=key_pos, dist=key_dist, vdist=v_dist, hdist=h_dist), 
                        x=graph.axis.linear(title=ind_title, min=min_ind, max=max_ind), 
                        y=graph.axis.linear(title=dep_title, min=min_dep, max=max_dep))
    
    # Create line styles that have predetermined color order for each pair in series.  
    # All Data Set 1, data is plotted with thick solid lines while Data Set 2, data is thin solid lines.
    # Colors are from http://pyx.sourceforge.net/manual/colorname.html
    d1PlotStyle = graph.style.line(lineattrs=[attr.changelist([color.cmyk.Black, color.cmyk.Red, color.cmyk.Green, color.cmyk.Blue]), style.linestyle.solid, style.linewidth(0.06*unit.w_cm)])
    d2PlotStyle = graph.style.line(lineattrs=[attr.changelist([color.cmyk.Grey, color.cmyk.Red, color.cmyk.Green, color.cmyk.Blue]), style.linestyle.solid, style.linewidth(0.03*unit.w_cm)])
    
    #Loop strcuture to process compound colum names in d line.
    if diagnostic_level >= 3:
        print "   <3> Data Set 1, Data:",d1_data
    if len(d1_data) > 1 : # or len(d2_data) > 1:
        #Set plot legend key text.
        d1_key_list = eval(plot_data['d1_Key'])
        d2_key_list = eval(plot_data['d2_Key'])
        d1_plot_counter = 0
        d2_plot_counter = 0
        
        # Loop through and plot Data Set 1, data
        for d1_data_item in d1_data:
            g.plot(graph.data.points(d1_data_item, title=d1_key_list[d1_plot_counter], x=1, y=2),
                  [d1PlotStyle])
            d1_plot_counter = d1_plot_counter + 1
            
        # Loop through and plot Data Set 1, data
        for d2_data_item in d2_data:
            g.plot(graph.data.points(d2_data_item, title=d2_key_list[d2_plot_counter], x=1, y=2),
                  [d2PlotStyle])
            d2_plot_counter = d2_plot_counter + 1
    else:
        #Set plot legend key text.
        d1_key = plot_data['d1_Key']
        d2_key = plot_data['d2_Key']
        
        if diagnostic_level >= 3:
            print "   <3> Data Set 1, Data to Plot:", d1_data[0]
            print "   <3> Data Set 2, Data to Plot:", d2_data[0]
        
        # Plot Data Set 1, data
        g.plot(graph.data.points(d1_data[0], title=d1_key, x=1, y=2),
            [d1PlotStyle])
        # Plot Predicted/Data Set 2, data
        g.plot(graph.data.points(d2_data[0], title=d2_key, x=1, y=2),
            [d2PlotStyle])
    
    # Now plot the Title text, alignment based on Title_Position setting.
    # "tl" or "tc" or "tr" or "ml" or "mc" or "mr" or "bl" or "bc" or "br"
    if title_position == 'tl':
        g.text(0.2, g.height - 0.2, plot_title, [text.halign.left,  text.valign.top, text.size.normalsize])
    elif title_position == 'tc':
        g.text(g.width/2, g.height - 0.2, plot_title, [text.halign.center, text.valign.top, text.size.normalsize])
    elif title_position == 'tr':
        g.text(g.width-0.2, g.height - 0.2, plot_title, [text.halign.right, text.valign.top, text.size.normalsize])
    elif title_position == 'ml':
        g.text(0.2, g.height/2, plot_title, [text.halign.left, text.valign.middle, text.size.normalsize])
    elif title_position == 'mc':
        g.text(g.width/2, g.height/2, plot_title, [text.halign.center, text.valign.middle, text.size.normalsize])
    elif title_position == 'mr':
        g.text(g.width-0.2, g.height/2, plot_title, [text.halign.right, text.valign.middle, text.size.normalsize])
    elif title_position == 'bl':
        g.text(0.2, 0.2, plot_title, [text.halign.left, text.valign.bottom, text.size.normalsize])
    elif title_position == 'bc':
        g.text(g.width/2, 0.2, plot_title, [text.halign.center, text.valign.bottom, text.size.normalsize])
    elif title_position == 'br':
        g.text(g.width-0.2, 0.2, plot_title, [text.halign.right, text.valign.bottom, text.size.normalsize])
    else:
        print "A title location was not specified.\nUsing the default top left position."
        g.text(0.2, g.height - 0.2, plot_title, [text.halign.left, text.valign.top, text.size.normalsize])
        
    # Write the output
    plot_file_path = output_directory+plot_file_name
    g.writePDFfile(plot_file_path)
    if diagnostic_level >= 2:
        print "\n*** Comparison Plot to: ***\n", plot_file_path+".PDF"

def scatter_plot(group_info,scatter_info,data_set):
    #data_set is a dictionary keyed by quantity, containing lists of groups and Independent Data and Dep data points.
    
    for quantity_number in scatter_info:
        #print "Dataset for quantity number "+str(quantity_number)+": ", data_set[quantity_number]
        
        if data_set[quantity_number] == []:
            if diagnostic_level >= 2:
                print "No Scatter Plot Data in Quantity "+str(quantity_number)+" Dataset.\n"
        else:
            if diagnostic_level >= 2:
                print "Scatter Plot Data for Quantity "+str(quantity_number)+" Dataset."
            # Set variables for Plot extracted from the first group of lines in config file starting with 'q'.
            
            # Variables for plot.
            plot_title = scatter_info[int(quantity_number)]['Scatter_Plot_Title']
            if diagnostic_level >= 2:
                print "Plot Title:", plot_title
            
            ind_title = scatter_info[int(quantity_number)]['Ind_Title']
            dep_title = scatter_info[int(quantity_number)]['Dep_Title']
            min_ind = float(scatter_info[int(quantity_number)]['Plot_Min'])
            max_ind = float(scatter_info[int(quantity_number)]['Plot_Max'])
            min_dep = float(scatter_info[int(quantity_number)]['Plot_Min'])
            max_dep = float(scatter_info[int(quantity_number)]['Plot_Max'])
            percent_error = float(scatter_info[int(quantity_number)]['%error'])
            title_position = scatter_info[int(quantity_number)]['Title_Position']
            plot_width = float(scatter_info[int(quantity_number)]['Plot_Width(cm)'])
            
            # Specify the position and line spacing of the plot key.
            key_pos  = scatter_info[int(quantity_number)]['Key_Position']
            key_dist = 0.1*unit.v_cm
            v_dist   = 0.7*unit.v_cm
            h_dist   = 0.4*unit.v_cm
            
            #Create filename from fields in input file record.
            plot_file_name = scatter_info[int(quantity_number)]['Plot_Filename']
            if diagnostic_level >= 3:
                print "   <3> Plot File Name:", plot_file_name
            
            # Determine the location for the key, alignment based on Key_Position setting.
            if key_pos == "tl" or "tc" or "tr" or "ml" or "mc" or "mr" or "bl" or "bc" or "br":
                pass
                if diagnostic_level >= 3:
                    print "   <3> Key Position =", key_pos
            else:
                if diagnostic_level >= 2:
                    print "The key position was not specified.\nUsing the default bottom right position."
                key_pos = "br"
            
            # Determine the location for the title, alignment based on Title_Position setting.
            if title_position == "tl" or "tc" or "tr" or "ml" or "mc" or "mr" or "bl" or "bc" or "br":
                ()
                if diagnostic_level >= 3:
                    print "   <3> Title Position =", title_position
            else:
                if diagnostic_level >= 2:
                    print "The Title position was not specified.\nUsing the default top left position."
                key_pos = "tl"
            
            #Begin Plotting
            
            # Initialize graph object
            g = graph.graphxy(width=plot_width, ratio=1/1, key=graph.key.key(pos=key_pos, dist=key_dist, vdist=v_dist, hdist=h_dist, textattrs=[text.size.footnotesize]), 
                                x=graph.axis.linear(title=ind_title, min=min_ind, max=max_ind), 
                                y=graph.axis.linear(title=dep_title, min=min_dep, max=max_dep))
            
            #Plot Midline and Error bounds lines.
            errorLineCenterPoints = [[min_ind,min_dep],[max_ind,max_dep]]
            if diagnostic_level >= 3:
                print "   <3> Error Line Center Points:", errorLineCenterPoints
            
            if min_ind < 0:       
                lower_bound = ((min_dep)+((min_dep)*(percent_error / 100)))
                errorLineLowerPoints = [[min_ind,lower_bound],[max_ind,max_dep]]
                upper_bound = ((min_dep)-((min_dep)*(percent_error/100)))
                errorLineUpperPoints = [[min_ind,upper_bound],[max_ind,max_dep]]
                if diagnostic_level >= 3:
                    print "   <3> Error Line Center Points:", errorLineCenterPoints
                    print "   <3> Lower Bound:", lower_bound
                    print "   <3> Lower Error Line Points:", errorLineLowerPoints
                    print "   <3> Upper Bound:", upper_bound
                    print "   <3> Upper Error Line Points:", errorLineUpperPoints
                                 
            else:
                lower_bound = max_dep - max_dep * percent_error / 100
                errorLineLowerPoints = [[min_ind,min_dep],[max_ind,lower_bound]]
                upper_bound = max_dep + max_dep * percent_error / 100.0
                errorLineUpperPoints = [[min_ind,min_dep],[max_ind,upper_bound]]
                if diagnostic_level >= 3:
                    print "   <3> Error Line Center Points:", errorLineCenterPoints
                    print "   <3> Lower Bound:", lower_bound
                    print "   <3> Lower Error Line Points:", errorLineLowerPoints
                    print "   <3> Upper Bound:", upper_bound
                    print "   <3> Upper Error Line Points:", errorLineUpperPoints
            
            g.plot(graph.data.points(errorLineCenterPoints, title=None, x=1, y=2),
                    [graph.style.line([style.linewidth.Thin, style.linestyle.solid])])
                    
            if percent_error == 0:
                if diagnostic_level >= 1:
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
            if diagnostic_level >= 3:
                print "   <3> Grouped Data List:", grouped_data_list
            
            if diagnostic_level >= 3:
                print "   <3> DataSet for Quantity "+str(quantity_number)+":", data_set[quantity_number]
            if len(data_set[quantity_number]) > 1:
                if diagnostic_level >= 3:
                    print "   <3> Grouped Scatter Data:", data_set[quantity_number]
                for arr_temp in grouped_data_list:
                    grouped_data_list[arr_temp] = []
                for data_set_item in data_set[quantity_number]:
                    if diagnostic_level >= 3:
                        print "   <3> Data Set Item:", data_set_item
                        print "   <3> Data for Group:"+data_set_item[0]+":", data_set_item[1]
                    grouped_data_list[int(data_set_item[0])].append(data_set_item[1])
                if diagnostic_level >= 3:
                    print "   <3> Grouped Data List:", grouped_data_list
                                
                group_counter = 0
                for j in grouped_data_list:
                    if diagnostic_level >= 3:
                        print "   <3> J =", j
                    if j != []:
                        if diagnostic_level >= 3:
                            print "   <3> Group Counter:", group_counter
                        
                        # Pull group symbol specifications from config file.
                        config_group_symbol = group_info[group_counter]["Symbol"]
                        if diagnostic_level >= 3:
                            print "   <3> Group Symbol:", config_group_symbol
                        group_symbol = "graph.style.symbol."+config_group_symbol
                        
                        config_group_symbol_color = group_info[group_counter]["Color"]
                        if diagnostic_level >= 3:
                            print "   <3> Group Symbol Color:", config_group_symbol_color
                        config_group_symbol_filled = group_info[group_counter]["Filled"]
                        if diagnostic_level >= 3:
                            print "   <3> Group Symbol Filled:", config_group_symbol_filled
                        
                        if config_group_symbol_filled == 'yes':
                            fillstyle = "deco.filled([color.cmyk."+config_group_symbol_color+"])"
                        else:
                            fillstyle = "deco.stroked([color.cmyk."+config_group_symbol_color+"])"  
                        
                        #Create temporary symbol style.
                        tempstyle = "graph.style.symbol("+group_symbol+", size=0.1*unit.v_cm, symbolattrs=["+fillstyle+"])"
                        
                        scatterpointstyle = eval(tempstyle)
                        
                        if diagnostic_level >= 3:
                            print "   <3> Group Title:", group_info[group_counter]["Group_Title"]
                        g.plot(graph.data.points(j, x=1, y=2, title=group_info[group_counter]["Group_Title"]), [scatterpointstyle])
                    else:
                        pass
                    group_counter = group_counter + 1
                    if diagnostic_level >= 3:
                        print "   <3> Group Counter:", group_counter
                    
            else:
                if diagnostic_level >= 1:
                    print "Non-Grouped Scatter Data:"
                scatter_plot_data = []
                scatter_plot_data.append(data_set[quantity_number][0][1])
                
            if diagnostic_level >= 3:
                print "   <3> Grouped Data List:", grouped_data_list
            
            # Now plot the Title text, alignment based on Title_Position setting.
            # "tl" or "tc" or "tr" or "ml" or "mc" or "mr" or "bl" or "bc" or "br"
            if title_position == 'tl':
                g.text(0.2, g.height - 0.2, plot_title, [text.halign.left,  text.valign.top, text.size.normalsize])
            elif title_position == 'tc':
                g.text(g.width/2, g.height - 0.2, plot_title, [text.halign.center, text.valign.top, text.size.normalsize])
            elif title_position == 'tr':
                g.text(g.width-0.2, g.height - 0.2, plot_title, [text.halign.right, text.valign.top, text.size.normalsize])
            elif title_position == 'ml':
                g.text(0.2, g.height/2, plot_title, [text.halign.left, text.valign.middle, text.size.normalsize])
            elif title_position == 'mc':
                g.text(g.width/2, g.height/2, plot_title, [text.halign.center, text.valign.middle, text.size.normalsize])
            elif title_position == 'mr':
                g.text(g.width-0.2, g.height/2, plot_title, [text.halign.right, text.valign.middle, text.size.normalsize])
            elif title_position == 'bl':
                g.text(0.2, 0.2, plot_title, [text.halign.left, text.valign.bottom, text.size.normalsize])
            elif title_position == 'bc':
                g.text(g.width/2, 0.2, plot_title, [text.halign.center, text.valign.bottom, text.size.normalsize])
            elif title_position == 'br':
                g.text(g.width-0.2, 0.2, plot_title, [text.halign.right, text.valign.bottom, text.size.normalsize])
            else:
                print "A title location was not specified.\nUsing the default top left position."
                g.text(0.2, g.height - 0.2, plot_title, [text.halign.left, text.valign.top, text.size.normalsize])
            
            #Make %error text on plot by error bars.
            # pos_percent_error = str(percent_error)+"%"
            # neg_percent_error = "-"+str(percent_error)+"%"
            # g.text(g.width - 0.4, g.height - 0.3, pos_percent_error, [text.halign.center, text.valign.middle, text.size.tiny])
            # g.text(g.width - 0.2, g.height - 0.4, neg_percent_error, [text.halign.center, text.valign.middle, text.size.tiny])
            
            # Write the output
            plot_file_path = output_directory+plot_file_name
            if diagnostic_level >= 3:
                print "   <3> Plot File Path:", plot_file_path
            g.writePDFfile(plot_file_path)
            if diagnostic_level >= 2:
                print "Scatter Plot to: \n", plot_file_path+".PDF\n"


### Start of Main Code
if diagnostic_level >= 1:
    print "**** READING CONFIGURATION FILE ****"

##Get information from config file.
group_quantity_data_dicts = extract_config_data(config_file_name)
if diagnostic_level >= 2:
    print "\nThere are "+str(len(group_quantity_data_dicts[0]))+" scatter data groups, (g lines)."
if diagnostic_level >= 3:
    print "   <3> Group Lines:\n", group_quantity_data_dicts[0]
if diagnostic_level >= 2:
    print "There are "+str(len(group_quantity_data_dicts[1]))+" scatter data quantities to plot, (q lines)."
if diagnostic_level >= 3:
    print "   <3> Quantity Lines:\n", group_quantity_data_dicts[1]
if diagnostic_level >= 2:
    print "There are "+str(len(group_quantity_data_dicts[2]))+" comparison data sets to plot, ("+data_line_char+" lines)."

## Create comparison plots
if diagnostic_level >= 1:
    print "**** CREATING COMPARISON DATA ****"
d_count = 1
for data_record in group_quantity_data_dicts[2]:
    # Each d line, data_record, may contain compound column names from the config file.
    if diagnostic_level >= 1:
        print "*** #"+str(d_count)+" of "+str(len(group_quantity_data_dicts[2]))+" comparison records. ***"
    # Extract relevant portions of comparison data as defined in config file.
    comp_data_to_plot = extract_comp_data(group_quantity_data_dicts[2][data_record])
    if diagnostic_level >= 3:
        print "   <3> Comparison Data to Plot:", comp_data_to_plot

    #Seperate d1 and d2 data lists.
    d1_plot_data = comp_data_to_plot[0]
    d2_plot_data = comp_data_to_plot[1]

    if diagnostic_level >= 3:
        print "   <3> Data Set 1, Plot Data:", d1_plot_data
        print "   <3> Data Set 2, Plot Data:", d2_plot_data

    # Create plot for data_record.
    if process_set == 1 or process_set == 2:
        comparison_plot(group_quantity_data_dicts[2][data_record],d1_plot_data,d2_plot_data)
    d_count = d_count + 1

## Create scatter plots
if process_set == 1 or process_set == 3:
    if diagnostic_level >= 1:
        print "**** CREATING SCATTER PLOTS ****"
    scatter_quantity = 1
    scatter_group = 1
    temp_scatter_data_list = []

    #Grouping Scatter Data by Quantity
    for scatter_plot_record in sorted(group_quantity_data_dicts[1]):
        if diagnostic_level >= 3:
            print "   <3> Scatter Plot:", scatter_plot_record
            print "   <3> Quantity:", group_quantity_data_dicts[1][scatter_plot_record]['Comparison_Quantities']
    
        for scatter_data_key in sorted(scatter_data_dict):
            split_key = split("~",scatter_data_key)
            if split_key[0] == str(scatter_plot_record) and scatter_data_dict[scatter_data_key] != []:
                temp_scatter_data_list.append([split_key[1], scatter_data_dict[scatter_data_key][:2]])
            else:
                pass

        combined_scatter_data[scatter_quantity] = temp_scatter_data_list
        temp_scatter_data_list = []
        scatter_quantity = scatter_quantity + 1

    # Plot Data
    if diagnostic_level >= 3:
        print "   <3> Data to Scatter Plot:", combined_scatter_data
    scatter_plot(group_quantity_data_dicts[0],group_quantity_data_dicts[1],combined_scatter_data)

## Write Summary Data to File.
#NRC Comparisons Output
# Output for each data set 
#*d1 Zero Val
#*d1 Peak Val
#*Peak Time Val
#*d2 Zero Val
#*d2 Peak Val
#*Peak Time Val
#*DeltaE
#*DeltaM
#*Rel Diff

if diagnostic_level >= 1:
    print "****  Processing finished, thank you for your patience.  ****"
