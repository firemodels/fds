import parser as pr
import calcs as calc
import plotter as plot

### Set Diagnostic Output Level
## Uncomment the level of diagnostics desired.
## 1 = Minimal, 2 = Normal, 3 = Maximum.
diagnostic_level = 2 # Input Value

print "**** Diagnostics Set at Level", diagnostic_level, "****\n"

### Set what plots to create: 
## BOTH (1), Comparison Only (2), Scatter Only (3)
process_set = 2 # Input Value

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

if diagnostic_level >= 1:
    print "**** Data Character Set to '"+data_line_char+"' ****\n"


## Validation Data Set (1), Verification Data Set (2), Examples Data Set (3), Trainier Data Set (4)
data_set = 5 # Input Value

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

## Examples Data 
if data_set == 3:
    data_directory = "../../Test_cases/"
    output_directory = "../../Manuals/FDS_5_User_Guide/FIGURES/"
    config_file_name = "Examples_Data_Config_File.csv"
    print "**** Processing Examples Data Set ****\n"

## Trainier Data
if data_set == 4:
    data_directory = "../../Training/"
    output_directory = "../../Manuals/FDS_SMV_Training_Guide/datafigures/"
    config_file_name = "Training_Examples_Data_Config_File.csv"
    print "**** Processing Training Data Set ****\n"
    
## Test Data
if data_set == 5:
    data_directory = "/Users/bwklein/Development/FDS-SMV/Validation/"
    output_directory = "./plots/"
    config_file_name = "validation_data_config"
    print "**** Processing Test Data Set ****\n"

### Set Global Variables to Zero or Null Sets.
config_files = ['groups','styles','quantities']

### Start of Main Code
if diagnostic_level >= 1:
    print "**** READING CONFIGURATION FILES ****"

##Get information from config files.
if diagnostic_level >= 2:
    print "Parsing Configs"
pr.parse_configs(config_files,diagnostic_level)

if diagnostic_level >= 2:
    print "Parsing Data Info"
pr.parse_data_info(config_file_name,data_line_char,diagnostic_level)

data_info = pr.read_pickle(config_file_name+'_object.pkl',diagnostic_level)
if diagnostic_level >= 1:
    print "There are",len(data_info),"data records to process."

#Start Processing records in data config file.
data_sets = pr.collect_data_sets(data_info,data_directory,diagnostic_level)

data_index_records = {} 
for key in data_sets:
    data_index_records[key] = {}
    data_index_records[key]['d1_index_set'] = pr.find_start_stop_index(data_sets[key][0][0],data_info[key]['d1_Start'],data_info[key]['d1_End'],data_info[key]['d1_Comp_Start'],data_info[key]['d1_Comp_End'],diagnostic_level)
    data_index_records[key]['d2_index_set'] = pr.find_start_stop_index(data_sets[key][1][0],data_info[key]['d2_Start'],data_info[key]['d2_End'],data_info[key]['d2_Comp_Start'],data_info[key]['d2_Comp_End'],diagnostic_level)
    #print data_index_records[key]['d1_index_set'],",",data_index_records[key]['d2_index_set']

# Make Comparison Plots.
if diagnostic_level >= 2:
    print "Begin Comparison Plots"

for key in data_info:
    #print key, data_info[key]['Plot_Filename']
    plot.comparison_plot(data_sets[key],data_info[key],data_index_records[key]['d1_index_set'],data_index_records[key]['d2_index_set'],output_directory,diagnostic_level)

if diagnostic_level >= 2:
    print "Finished Comparison Plots"

# Create Scatter Plot Data Object based on quantities object.
pr.make_scatter_dict(diagnostic_level)

# Find metric data for scatter plots and write data to scatter_data_dict object.
for key in data_info:
    if data_info[key]['Quantity'] == str(0):
        if diagnostic_level >= 1:
            print "Quantity set to 0, no scatter data."
            #This allows the d line Quantity value to be set to 0 when either d1 or d2 data is missing.
    else:
        d1_metric_start_index = data_index_records[key]['d1_index_set'][2]
        d1_metric_stop_index = data_index_records[key]['d1_index_set'][3]
        d2_metric_start_index = data_index_records[key]['d2_index_set'][2]
        d2_metric_stop_index = data_index_records[key]['d2_index_set'][3]
        
        if data_info[key]['Metric'] == 'max':
            for data_index in range(len(data_sets[key][0][1:])):
                max_results = calc.calc_max(data_sets[key][0][data_index+1][d1_metric_start_index:d1_metric_stop_index],data_sets[key][1][data_index+1][d2_metric_start_index:d2_metric_stop_index],data_info[key]['d1_Initial_Value'],data_info[key]['d2_Initial_Value'],diagnostic_level)
                pr.add_group_data_to_scatter_data_dict(data_info[key]['Quantity'],data_info[key]['Group'],max_results[:2],diagnostic_level)
        elif data_info[key]['Metric'] == 'min':
            for data_index in range(len(data_sets[key][0][1:])):
                min_results = calc.calc_min(data_sets[key][0][data_index+1][d1_metric_start_index:d1_metric_stop_index],data_sets[key][1][data_index+1][d2_metric_start_index:d2_metric_stop_index],data_info[key]['d1_Initial_Value'],data_info[key]['d2_Initial_Value'],diagnostic_level)
                pr.add_group_data_to_scatter_data_dict(data_info[key]['Quantity'],data_info[key]['Group'],min_results[:2],diagnostic_level)
        else:
            if diagnostic_level >= 1:
                print "!!! Metric is undefined in the input file. !!!"
            exit()

# Make Scatter Plots.
if diagnostic_level >= 1:
    print "Begin Scatter Plots"

scatter_data_dict = pr.read_pickle('scatter_data_dict_object.pkl',diagnostic_level)
for quantity_key in scatter_data_dict:
    if scatter_data_dict[quantity_key] != {}:
        plot.scatter_plot(quantity_key,scatter_data_dict[quantity_key],output_directory,diagnostic_level)

if diagnostic_level >= 1:
    print "Finished Scatter Plots"


## Report Completion.
if diagnostic_level >= 1:
    print "****  Processing finished, thank you for your patience.  ****"