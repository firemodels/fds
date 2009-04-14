import os
import modules.pyrograph_parser as prsr
import modules.pyrograph_calcs as calc
import modules.pyrograph_plotter as plot
import modules.pyrograph_progress as prgrs

import time
start_time = time.time()

print "\n***** PyroGraph Data Processing Utility *****\n"

#Process User Configuration file.
user_config = {}
 
file_name = "pyrograph_config.txt"
config_file= open(file_name)
for line in config_file:
    line = line.strip()
    if line and line[0] is not "#" and line[-1] is not "=":
        var,val = line.rsplit("=",1)
        user_config[var.strip()] = val.strip()
config_file.close()

config_files = ['pyrograph_groups','pyrograph_styles','pyrograph_quantities']

diagnostic_level = int(user_config['diagnostic_level'])
print "**** Diagnostics Set at Level", diagnostic_level

process_set = int(user_config['process_set'])
if diagnostic_level >= 1:
    if process_set == 1:
        print "**** Plotting both Comparison and Scatter Data"
    if process_set == 2:
        print "**** Plotting Only Comparison Data"
    if process_set == 3:
        print "**** Plotting Only Scatter Data"

data_line_char = user_config['data_line_char']

if diagnostic_level >= 1:
    print "**** Data Character Set to '"+data_line_char+"'"

data_set = int(user_config['data_set'])

## Validation Data
if data_set == 1:
    data_directory = "../../Validation/"
    output_directory = "../../Manuals/FDS_5_Validation_Guide/FIGURES/"
    config_file_name = "validation_data_config"
    print "\n** Processing Validation Data Set\n"

## Verification Data
if data_set == 2:
    data_directory = "../../Verification/"
    output_directory = "../../Manuals/"
    config_file_name = "verification_data_config"
    print "\n** Processing Verification Data Set\n"

## Trainier Data
if data_set == 4:
    data_directory = "../../Training/"
    output_directory = "../../Manuals/FDS_SMV_Training_Guide/datafigures/"
    config_file_name = "Training_Examples_Data_Config_File.csv"
    print "\n** Processing Training Data Set\n"


### Start of Main Code

## Get information from pyrograph configuration files.
if diagnostic_level >= 2:
    print "Parsing Configs"
prsr.parse_configs(config_files,diagnostic_level)

if diagnostic_level >= 2:
    print "Parsing Data Info"
prsr.parse_data_info(config_file_name,data_line_char,diagnostic_level)

data_info = prsr.read_pickle(config_file_name+'_object.pkl',diagnostic_level)
if diagnostic_level >= 2:
    print "There are",len(data_info),"data records to process."

#Start Processing records in data config file.
data_sets = prsr.collect_data_sets(data_info,data_directory,diagnostic_level)
styles = prsr.read_pickle("pyrograph_styles_object.pkl",diagnostic_level)
quantities = prsr.read_pickle("pyrograph_quantities_object.pkl",diagnostic_level)
groups = prsr.read_pickle("pyrograph_groups_object.pkl",diagnostic_level)

data_index_records = {}

for key in data_sets:
    data_index_records[key] = {}
    data_index_records[key]['d1_index_set'] = prsr.find_start_stop_index(data_sets[key][0][0], \
    data_info[key]['d1_Start'],data_info[key]['d1_End'],data_info[key]['d1_Comp_Start'], \
    data_info[key]['d1_Comp_End'],data_info[key]['Scale_Ind'],diagnostic_level)
    
    data_index_records[key]['d2_index_set'] = prsr.find_start_stop_index(data_sets[key][1][0], \
    data_info[key]['d2_Start'],data_info[key]['d2_End'],data_info[key]['d2_Comp_Start'], \
    data_info[key]['d2_Comp_End'],data_info[key]['Scale_Ind'],diagnostic_level)

# Comparison Plotting...
if process_set == 1 or process_set == 2:
    if diagnostic_level >= 2:
        print "Begin Comparison Plots"
    
    if diagnostic_level < 3:
        total = len(data_sets)
        p = prgrs.ProgressMeter(total=total, unit='Comparison Plots', rate_refresh=0.25)
    
    for key in data_info:
        plot.comparison_plot(data_sets[key],data_info[key],data_index_records[key]['d1_index_set'], \
        data_index_records[key]['d2_index_set'],styles,output_directory,diagnostic_level)
        
        if diagnostic_level < 3:
            p.update(1)
        
    if diagnostic_level >= 2:
        print "Finished Comparison Plots"

# Scatter Plotting...
if process_set == 1 or process_set == 3:
    # Create Scatter Plot Data Object based on quantities object.
    prsr.make_scatter_dict(config_files[2],diagnostic_level)

    # Find metric data for scatter plots and write data to scatter_data_dict object.
    for key in data_info:
        if data_info[key]['Quantity'] == str(0):
            if diagnostic_level >= 3:
                print "Quantity set to 0, no scatter data."
                #This allows the d line Quantity value to be set to 0 when either d1 or d2 data is missing.
        else:
            d1_metric_start_index = data_index_records[key]['d1_index_set'][2]
            d1_metric_stop_index = data_index_records[key]['d1_index_set'][3]
            d2_metric_start_index = data_index_records[key]['d2_index_set'][2]
            d2_metric_stop_index = data_index_records[key]['d2_index_set'][3]
            
            if diagnostic_level >= 3:
                print "Processing Data for:",data_info[key]['Plot_Filename']
                
            if data_info[key]['Metric'] == 'max':
                for data_index in range(len(data_sets[key][0][1:])):
                    max_results = calc.calc_max(data_sets[key][0][data_index+1][d1_metric_start_index:d1_metric_stop_index], \
                    data_sets[key][1][data_index+1][d2_metric_start_index:d2_metric_stop_index],data_info[key]['d1_Initial_Value'], \
                    data_info[key]['d2_Initial_Value'],diagnostic_level)
                    
                    prsr.add_group_data_to_scatter_data_dict(data_info[key]['Quantity'],data_info[key]['Group'],max_results[:2],diagnostic_level)
            elif data_info[key]['Metric'] == 'min':
                for data_index in range(len(data_sets[key][0][1:])):
                    min_results = calc.calc_min(data_sets[key][0][data_index+1][d1_metric_start_index:d1_metric_stop_index], \
                    data_sets[key][1][data_index+1][d2_metric_start_index:d2_metric_stop_index],data_info[key]['d1_Initial_Value'], \
                    data_info[key]['d2_Initial_Value'],diagnostic_level)
                    
                    prsr.add_group_data_to_scatter_data_dict(data_info[key]['Quantity'],data_info[key]['Group'],min_results[:2],diagnostic_level)
            else:
                if diagnostic_level >= 1:
                    print "!!! Metric is undefined in the input file. !!!"
                exit()
                
    scatter_data_dict = prsr.read_pickle('scatter_data_dict_object.pkl',diagnostic_level)
    dict_length = 0
    
    for key in scatter_data_dict:
        dict_length = dict_length + len(scatter_data_dict[key])

    if dict_length > 0:
        # Make Scatter Plots.
        if diagnostic_level >= 1:
            print "Begin Scatter Plots"
            
        count_no_data = 0
        for quantity_key in scatter_data_dict:
            if scatter_data_dict[quantity_key] == {}:
                count_no_data += 1
        if diagnostic_level < 3:
            total = len(scatter_data_dict.keys())-count_no_data
            p = prgrs.ProgressMeter(total=total, unit='Scatter Plots', rate_refresh=0.25)
    
        for quantity_key in scatter_data_dict:
            if scatter_data_dict[quantity_key] != {}:
                plot.scatter_plot(quantity_key,scatter_data_dict[quantity_key],quantities,groups,styles,output_directory,diagnostic_level)
                if diagnostic_level < 3:
                    p.update(1)
            
        if diagnostic_level >= 1:
            print "Finished Scatter Plots"
    else:
        print "No Scatter Plots to Create, all Quantities set to 0.\nConsider setting process_set=2 in pyrograph_config.txt"

# Clean up temp files
for file_name in config_files:
    os.remove(file_name+'_object.pkl')
if process_set == 1 or process_set == 3:
    os.remove('scatter_data_dict_object.pkl')
os.remove(config_file_name+'_object.pkl')

## Report Completion.
if diagnostic_level >= 1:
    print "****  Processing finished"
    run_time = time.time() - start_time
    print '      Run time: %f seconds' % run_time
