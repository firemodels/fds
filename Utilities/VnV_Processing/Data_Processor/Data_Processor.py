import parser as pr
import calcs as calc
import plotter as plot

### Set Diagnostic Output Level
## Uncomment the level of diagnostics desired.
## 1 = Minimal, 2 = Normal, 3 = Maximum.
diagnostic_level = 1 # Input Value

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
if diagnostic_level >= 2:
    print "**** READING CONFIGURATION FILES ****"

##Get information from config files.
print "Parsing Configs"
pr.parse_configs(config_files,diagnostic_level)
print "Parsing Data Info"
pr.parse_data_info(config_file_name,data_line_char,diagnostic_level)
print "Finished Parsing Files"

data_info = pr.read_pickle(config_file_name+'_object.pkl',diagnostic_level)
print "There are",len(data_info),"data records to process."

#Start Processing records in data config file.
data_sets = pr.collect_data_sets(data_info,data_directory,diagnostic_level)

data_index_records = {}
for key in data_sets:
    data_index_records[key] = {}
    #print "Find Start/Stop for d1 in", data_info[key]['d1_Filename']
    data_index_records[key]['d1_index_set'] = pr.find_start_stop_index(data_sets[key][0][0],data_info[key]['d1_Start'],data_info[key]['d1_End'],data_info[key]['d1_Comp_Start'],data_info[key]['d1_Comp_End'],data_info[key]['Scale_Ind'],diagnostic_level)
    #print "Find Start/Stop for d2 in", data_info[key]['d2_Filename']
    data_index_records[key]['d2_index_set'] = pr.find_start_stop_index(data_sets[key][1][0],data_info[key]['d2_Start'],data_info[key]['d2_End'],data_info[key]['d2_Comp_Start'],data_info[key]['d2_Comp_End'],data_info[key]['Scale_Ind'],diagnostic_level)
    #print data_index_records[key]['d1_index_set'],",",data_index_records[key]['d2_index_set']
print "Number of start/stop records",len(data_index_records.keys())

# Create Scatter Plot Data Object
pr.make_scatter_dict(diagnostic_level)

#This allows the d line Quantity value to be set to 0 when either d1 or d2 data is missing.

for key in data_info:
    if data_info[key]['Quantity'] == str(0):
        if diagnostic_level >= 2:
            print "Quantity set to 0, no comparison made."
    else:
        d1_metric_start_index = data_index_records[key]['d1_index_set'][2]
        d1_metric_stop_index = data_index_records[key]['d1_index_set'][3]
        d2_metric_start_index = data_index_records[key]['d2_index_set'][2]
        d2_metric_stop_index = data_index_records[key]['d2_index_set'][3]
        
        #print "Indices:", d1_metric_start_index, d1_metric_stop_index, d2_metric_start_index, d2_metric_stop_index
        
        if data_info[key]['Metric'] == 'max':
            #print "Max Calc..."
            for data_index in range(len(data_sets[key][0][1:])):
                max_results = calc.calc_max(data_sets[key][0][data_index+1][d1_metric_start_index:d1_metric_stop_index],data_sets[key][1][data_index+1][d2_metric_start_index:d2_metric_stop_index],data_info[key]['d1_Initial_Value'],data_info[key]['d2_Initial_Value'],diagnostic_level)
                #Add Results to Scatter Data Dictionary.
                #print "Data Index:", data_index
                #print "Max Results:", max_results[:2]
                pr.add_group_data_to_scatter_data_dict(data_info[key]['Quantity'],data_info[key]['Group'],max_results[:2],diagnostic_level)
        elif data_info[key]['Metric'] == 'min':
            #print "Min Calc..."
            for data_index in range(len(data_sets[key][0][1:])):
                min_results = calc.calc_min(data_sets[key][0][data_index+1][d1_metric_start_index:d1_metric_stop_index],data_sets[key][1][data_index+1][d2_metric_start_index:d2_metric_stop_index],data_info[key]['d1_Initial_Value'],data_info[key]['d2_Initial_Value'],diagnostic_level)
                #Add Results to Scatter Data Dictionary.
                #print "Data Index:", data_index
                #print "Min Results:", min_results[:2]
                pr.add_group_data_to_scatter_data_dict(data_info[key]['Quantity'],data_info[key]['Group'],min_results[:2],diagnostic_level)
        else:
            print "!!! Metric is undefined in the input file. !!!"
            exit()

scatter_data_dict = pr.read_pickle('scatter_data_dict_object.pkl',diagnostic_level)
for quantity_key in scatter_data_dict:
    print quantity_key
    if scatter_data_dict[quantity_key] != {}:
        plot.scatter_plot(quantity_key,scatter_data_dict[quantity_key],output_directory,diagnostic_level)
    #scatter_data_dict[quantity_key]



## Call stats functions, which write Summary of results to csv file.
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


## Call Plotter for Comparison Data.
if diagnostic_level >= 1:
    print "**** CREATING COMPARISON DATA ****"

d_count = 1
# 
# for data_record in group_quantity_data_dicts[2]:
#     # Each d line, data_record, may contain compound column names from the config file.
#     if diagnostic_level >= 1:
#         print "*** #"+str(d_count)+" of "+str(len(group_quantity_data_dicts[2]))+" comparison records. ***"
#     # Extract relevant portions of comparison data as defined in config file.
#     comp_data_to_plot = extract_comp_data(group_quantity_data_dicts[2][data_record])
#     if diagnostic_level >= 3:
#         print "Comparison Data to Plot:", comp_data_to_plot
# 
#     #Seperate d1 and d2 data lists.
#     d1_plot_data = comp_data_to_plot[0]
#     d2_plot_data = comp_data_to_plot[1]
# 
#     if diagnostic_level >= 3:
#         print "Data Set 1, Plot Data:", d1_plot_data
#         print "Data Set 2, Plot Data:", d2_plot_data
# 
#     # Create plot for data_record.
#     if process_set == 1 or process_set == 2:
#         comparison_plot(group_quantity_data_dicts[2][data_record],d1_plot_data,d2_plot_data)
#     d_count = d_count + 1
# 
# ## Call Plotter for Scatter Data.
# if process_set == 1 or process_set == 3:
#     if diagnostic_level >= 1:
#         print "**** CREATING SCATTER PLOTS ****"
#     scatter_quantity = 1
#     scatter_group = 1
#     temp_scatter_data_list = []
# 
#     #Grouping Scatter Data by Quantity
#     for scatter_plot_record in sorted(group_quantity_data_dicts[1]):
#         if diagnostic_level >= 3:
#             print "Scatter Plot:", scatter_plot_record
#             print "Quantity:", group_quantity_data_dicts[1][scatter_plot_record]['Comparison_Quantities']
#     
#         for scatter_data_key in sorted(scatter_data_dict):
#             split_key = split("~",scatter_data_key)
#             if split_key[0] == str(scatter_plot_record) and scatter_data_dict[scatter_data_key] != []:
#                 temp_scatter_data_list.append([split_key[1], scatter_data_dict[scatter_data_key][:2]])
#             else:
#                 pass
# 
#         combined_scatter_data[scatter_quantity] = temp_scatter_data_list
#         temp_scatter_data_list = []
#         scatter_quantity = scatter_quantity + 1
# 
#     # Plot Data
#     if diagnostic_level >= 3:
#         print "Data to Scatter Plot:", combined_scatter_data
#     scatter_plot(group_quantity_data_dicts[0],group_quantity_data_dicts[1],combined_scatter_data)
# 

## Report Completion.
if diagnostic_level >= 1:
    print "****  Processing finished, thank you for your patience.  ****"