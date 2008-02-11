import csv
from re import *
from pyx import *
from decimal import Decimal

### Set Global Variables

data_directory = "../../Validation/"
config_file_name = "Validation_Data_Config_File.csv"

scatter_data_dict = {}
combined_scatter_data = {}

### Define Functions

def extract_config_data(config_file):
    # Collect Data from Config File
    quantities_dict = {}
    data_dict = {}
    
    fh = file(config_file, 'U')
    data_array = csv.reader(fh)
    
    config_lists = [list(sublist) for sublist in data_array]
    #print config_test_lists
    print str(len(config_lists)) + " lines read from Configuration file\n"

    for each_list in config_lists:
        #print "First field is "+each_list[0]
        if each_list[0] == 'q':
            quantities_dict[each_list[1]] = each_list[2:]
        elif each_list[0] == 'd':
            data_dict[each_list[2]+"-"+each_list[4]] = each_list[1:]
        else:
            print """No information was found for quantities or data./nDouble check config file name."""
    
    # Remove header entries from Dictionaries, maybe pop these out to lists for use later.
    quant_header = quantities_dict.pop('Index')
    data_header = data_dict.pop('Dataname-Exp_Col_Name')
    
    #Convert quantities keys from string to integer.
    tempdict = dict((int(k), v) for k,v in quantities_dict.iteritems())
    quantities_dict = tempdict
    
    #Diagnostics:
    #print sorted(quantities_dict)
    #print quant_header
    #print data_header
    
    # Return a single list object containing the dictionaries.
    return [quant_header,quantities_dict,data_header,data_dict]

def extract_comp_data(comp_file_info):
    ## Read in from config file and Process data from Source .csv files.
    
    #print comp_file_info
    
    exp_data_dict = {}
    mod_data_dict = {}
    
    #List of variables from configuration file
    
    exp_data_filename = comp_file_info[2] #String of filename
    exp_column_name = comp_file_info[3] #Experimental Data Column Name
    mod_data_filename = comp_file_info[10] #String of filename
    mod_column_name = comp_file_info[11] #Experimental Data Column Name
    
    scatter_data_label = comp_file_info[0]+"-"+comp_file_info[1]+"-"+comp_file_info[3]
    
    exp_start_time_data_val = comp_file_info[5] #String in minutes to start exp plot data
    exp_stop_time_data_val = comp_file_info[6]  #String in minutes to stop exp plot data
    exp_start_time_comp_val = comp_file_info[7] #String in minutes to start exp compare data
    exp_stop_time_comp_val = comp_file_info[8]  #String in minutes to start exp compare data
    exp_initial_value = comp_file_info[9]       #Initial Value for Quantity
    
    mod_start_time_data_val = comp_file_info[13] #String in minutes to start mod plot data
    mod_stop_time_data_val = comp_file_info[14]  #String in minutes to stop mod plot data
    mod_start_time_comp_val = comp_file_info[15] #String in minutes to start mod compare data
    mod_stop_time_comp_val = comp_file_info[16]  #String in minutes to start mod compare data
    mod_initial_value = comp_file_info[17]       #Initial Value for Quantity
    
    min_max = comp_file_info[18] #String indicating if min or max value is required.

    def find_start_stop_index(data_dict,col_name,start_time_data,stop_time_data,start_time_comp,stop_time_comp):
        #This function is used to find index numbers for start and stop points in plotting and min-max values.
        
        ##Echo Input Values Diagnostic.
        #print start_time_data
        #print stop_time_data
        #print start_time_comp
        #print stop_time_comp
        
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
            if time_value2 < (float(stop_time_data)*60):
                rowcounter2 += 1
            else:
                row_number2 = (rowcounter2 - 1)
                #print "Set #2"
                #print "Time Ends at row #: "+str(row_number2)
                #print "With a value of: "+str(data_dict[col_name][row_number2])
                time_end_index = row_number2
                break
                
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

    exp_file_object = open(data_directory+exp_data_filename, "rb")
    mod_file_object = open(data_directory+mod_data_filename, "rb")
    
    ## Start File Processing
    
    #Read in experimental data and flip lists from rows to columns.
    exp_data_cols = zip(*csv.reader(exp_file_object))
    #Convert tuples to lists.
    exp_data_list = [list(sublist) for sublist in exp_data_cols]
    #Pull the time Column name out and strip whitespace from ends of string.
    exp_time_col_name = (exp_data_list[0][0]).strip()
    
    #Build Experimental Data Dictionary
    for exp_list in exp_data_list:
        exp_data_dict[(exp_list[0]).strip()] = map(float, exp_list[1:])
    
    #Read in model data and flip lists from rows to columns.
    mod_data_cols = zip(*csv.reader(mod_file_object))
    #Convert tuples to lists.
    mod_data_list = [list(sublist) for sublist in mod_data_cols]
    #Pull the time Column name out and strip whitespace from ends of string.
    mod_time_col_name = (mod_data_list[0][0]).strip()
    
    #Build Prediction/Model Data Dictionary
    for mod_list in mod_data_list:
        mod_data_dict[(mod_list[0]).strip()] = map(float, mod_list[1:])

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

def compute_difference(m_p,m_o,e_p,e_o):
    # Equation to compute relative difference of model and experimental data.
    M_p = float(m_p) #Peak value of Model Prediction to float
    M_o = float(m_o) #Original value of Model Prediction to float
    E_p = float(e_p) #Peak value of Experimental Measurement to float
    E_o = float(e_o) #Original value of Experimental Measurement to float
    e = (((M_p-M_o)-(E_p-E_o))/(E_p-E_o))
    return e


## Plot Validation Data to PDF files
## Two kinds of plots... Comparison and Scatter

def comparison_plot(plot_data,exp_data,mod_data):
    #plot_data is the items from the 'd' rows of the config file.
    
    # Variables for plot.
    plot_title = plot_data[19]
    x_title = plot_data[20]
    y_title = plot_data[21]   
    min_x = float(plot_data[22])
    max_x = float(plot_data[23])
    min_y = float(plot_data[24])
    max_y = float(plot_data[25])
    title_quadrant = int(plot_data[26])
    key_quadrant = int(plot_data[27])
    plot_width = int(plot_data[28])
    
    #Set plot legend key text.
    exp_key = plot_data[4]
    mod_key = plot_data[12] 
    
    #Create filename from fields in input file record.
    plot_file_name = plot_data[1]+"-"+plot_data[3]+".pdf"
    print plot_file_name
    
    # Determine the location for the key, alignment based on key_quadrant setting.
    if key_quadrant == 1:
        key_pos = "tl"
        key_dist = 0.00
    elif key_quadrant == 2:
        key_pos = "tr"
        key_dist = 0.00
    elif key_quadrant == 3:
        key_pos = "bl"
        key_dist = 0.00
    elif key_quadrant == 4:
        key_pos = "br"
        key_dist = 0.00
    else:
        print "A quadrant for the key location was not specified./nUsing the default top left quadrant."
        key_pos = "tl"
        key_dist = 0.00
    
    #Begin Plotting
    #print exp_data
    #print mod_data
    # Initialize graph object
    g = graph.graphxy(width=plot_width, ratio=4./3, key=graph.key.key(pos=key_pos, dist=key_dist), 
                        x=graph.axis.linear(title=x_title, min=min_x, max=max_x), 
                        y=graph.axis.linear(title=y_title, min=min_y, max=max_y))

    #One line at a time added to plot from each data set.
    # Plot Experimental data
    g.plot(graph.data.points(exp_data, title=exp_key, x=1, y=2),
        [graph.style.line([style.linewidth.Thin, style.linestyle.solid])])
    # Plot Predicted/Model data
    g.plot(graph.data.points(mod_data, title=mod_key, x=1, y=2),
        [graph.style.line([style.linewidth.Thin, style.linestyle.dashed])])
    
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
    g.writePDFfile(plot_file_name)


def scatter_plot(plot_info,data_set):
    #plot_info is details about the overall plot layout and text.
    #data_set is a dictionary keyed by test name containing lists of X and Y data points.
    print plot_info
    print data_set
    
    for quantity_number in plot_info:
        # Set variables for Plot extracted from the first group of lines in config file starting with 'q'.
        
        # Variables for plot.
        plot_title = plot_info[int(quantity_number)][1]
        print plot_title
        x_title = plot_info[int(quantity_number)][2]
        y_title = plot_info[int(quantity_number)][3]
        min_x = float(plot_info[int(quantity_number)][4])
        print min_x
        max_x = float(plot_info[int(quantity_number)][5])
        print max_x
        min_y = float(plot_info[int(quantity_number)][4])
        max_y = float(plot_info[int(quantity_number)][5])
        percent_error = int(plot_info[int(quantity_number)][6])
        title_quadrant = int(plot_info[int(quantity_number)][7])
        key_quadrant = int(plot_info[int(quantity_number)][8])
        plot_width = int(plot_info[int(quantity_number)][9])
        
        #Create filename from fields in input file record.
        plot_file_name = plot_info[int(quantity_number)][0]+".pdf"
        print plot_file_name
        
        # Determine the location for the key, alignment based on key_quadrant setting.
        if key_quadrant == 1:
            key_pos = "tl"
            key_dist = 0.00
        elif key_quadrant == 2:
            key_pos = "tr"
            key_dist = 0.00
        elif key_quadrant == 3:
            key_pos = "bl"
            key_dist = 0.00
        elif key_quadrant == 4:
            key_pos = "br"
            key_dist = 0.00
        else:
            print "A quadrant for the key location was not specified./nUsing the default top left quadrant."
            key_pos = "tl"
            key_dist = 0.00
        
        #Begin Plotting
        #print exp_data
        #print mod_data
        # Initialize graph object
        g = graph.graphxy(width=plot_width, ratio=1/1, key=graph.key.key(pos=key_pos, dist=key_dist), 
                            x=graph.axis.linear(title=x_title, min=min_x, max=max_x), 
                            y=graph.axis.linear(title=y_title, min=min_y, max=max_y))
        
        #Plot Midline and Error bounds lines.
        errorLineCenterPoints = [[min_x,min_y],[max_x,max_y]]
        print errorLineCenterPoints
        lower_bound = max_y - max_y * percent_error / 100
        print lower_bound
        errorLineLowerPoints = [[min_x,min_y],[max_x,lower_bound]]
        print errorLineLowerPoints
        upper_bound = max_y + max_y * percent_error / 100.0
        print upper_bound
        errorLineUpperPoints = [[min_x,min_y],[max_x,upper_bound]]
        print errorLineUpperPoints
        
        g.plot(graph.data.points(errorLineCenterPoints, title=None, x=1, y=2),
                [graph.style.line([style.linewidth.Thin, style.linestyle.solid])])
                
        g.plot(graph.data.points(errorLineLowerPoints, title=None, x=1, y=2),
                [graph.style.line([style.linewidth.Thin, style.linestyle.dashed])])
                
        g.plot(graph.data.points(errorLineUpperPoints, title=None, x=1, y=2),
                [graph.style.line([style.linewidth.Thin, style.linestyle.dashed])])
    
        #One line at a time added to plot from each data set.
        # Iterate over items in data dictionary key for keys that are not [].
        for data in data_set[quantity_number]:
        # Plot Experimental data
            print data[0]
            print data[1]
            #g.plot(graph.data.points(data[1], title=data[0], x=1, y=2))
            #    [graph.style.symbol()])
  
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
        g.writePDFfile(plot_file_name)
    

### Start of Main Code

#Get information from config file.
config_and_data_dicts = extract_config_data(config_file_name)
print "There are "+str(len(config_and_data_dicts[1]))+" quantities to process."
print "There are "+str(len(config_and_data_dicts[3]))+" data sets to plot."

## Create comparison plots
for data_record in config_and_data_dicts[3]:
    # Extract relevant portions of comparison data as defined in config file.
    comp_data_to_plot = extract_comp_data(config_and_data_dicts[3][data_record])
    exp_plot_data = comp_data_to_plot[0]
    mod_plot_data = comp_data_to_plot[1]
    
    # Create plots
    #comparison_plot(config_and_data_dicts[3][data_record],exp_plot_data,mod_plot_data)


## Create scatter plots
scatter_quantity = 1
temp_scatter_data_list = []

for scatter_plot_record in config_and_data_dicts[1]:
    #print "s record", scatter_plot_record
    #print "Quantity:", config_and_data_dicts[1][scatter_plot_record][0]
    for scatter_data_key in sorted(scatter_data_dict):
    	#print scatter_quantity
        split_key = split("-",scatter_data_key)
        if scatter_data_key[0] == str(scatter_plot_record) and scatter_plot_record != []:
            #print "Test name:", split_key[1]
            #print scatter_data_dict[scatter_data_key][:2]
            temp_scatter_data_list.append([split_key[1], scatter_data_dict[scatter_data_key][:2]])
            
    combined_scatter_data[scatter_quantity] = temp_scatter_data_list
    temp_scatter_data_list = []
    scatter_quantity = scatter_quantity + 1
        
    #scatter_plot_info = config_and_data_dicts[1][scatter_plot_record]
    #print scatter_plot_info	print scatter_data_key, scatter_data_dict[scatter_data_key]

#print combined_scatter_data
    
# Create plots
scatter_plot(config_and_data_dicts[1],combined_scatter_data)
    
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