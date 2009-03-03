"""
parser.py

This parser module processes the config files it is passed and 
then writes the resulting three data sets to pickle objects.

1) Comparison Data for Plotting
2) Data for Analytical Processing
3) Scatter Data for Plotting

"""

import sys
import os
import csv
import cPickle as cP

#diagnostic_level = 1
config_files = ['groups','styles','quantities']
data_config_file = 'validation_data_config'
scatter_data_dict = {}

def read_config(config_file,diagnostic_level):
    #Create File Object
    try:
        fh = file(config_file, 'U')
    except:
        print"!!! The Config File "+config_file+" does not exist or the path defined in the script is incorrect. !!!"
        exit()
    
    #Read file with DictReader from csv module.
    reader = csv.DictReader(fh, delimiter=",")
    data = {}
    while True:
        try:
            rdr = reader.next()
            ID_num = int(rdr['ID'])
            data[ID_num] = rdr
            # Remove redundant ID values from dictionaries
            data[ID_num].pop('ID')
        except StopIteration: break
    
    if diagnostic_level >= 3:
        print "Config Data:", data
    fh.close()
    # Return Configuration dictionary object.
    return data

def read_data_config(config_file,data_line_char,diagnostic_level):
    #Create File Object
    try:
        fh = file(config_file, 'U')
    except:
        print"!!! The Config File "+config_file+" does not exist or the path defined in the script is incorrect. !!!"
        exit()
    
    #Read file with DictReader from csv module.
    reader = csv.DictReader(fh, delimiter=",")
    data = {}
    record_count = 0
    while True:
        try:
            rdr = reader.next()
            # print rdr
            if rdr['switch_id'] == 'e':
                print "e character detected, stopping now."
                break
            elif rdr['switch_id'] == data_line_char:
                record_count += 1
                data_key_name = str(record_count)+"~"+rdr['Quantity'].strip()+"~"+rdr['Group'].strip()
                data[data_key_name] = rdr
                # Remove redundant ID values from dictionaries
                data[data_key_name].pop('switch_id')
            else:
                pass
                #print "Skipping Row."
        except StopIteration: break
    if diagnostic_level >= 3:
        print "Data Info:", data
    fh.close()
    # Return Configuration dictionary object.
    return data

def pickle_object(name,data_object,diagnostic_level):
    fh = file(name+'_object.pkl', 'w')
    cP.dump(data_object, fh, protocol=2)
    fh.close()
    if diagnostic_level >= 3:
        print "Wrote to "+name+"_object.pkl"

def read_pickle(pickle_file_name,diagnostic_level):
    fh = file(pickle_file_name, 'rb')
    data = cP.load(fh)
    fh.close()
    if diagnostic_level >= 3:
        print "Read information from "+pickle_file_name
    return data

def parse_configs(config_files,diagnostic_level):
    for config_name in config_files:
        config_data = read_config(config_name+'.csv',diagnostic_level)    
        pickle_object(config_name,config_data,diagnostic_level)
        if diagnostic_level >= 3:
            print "Parsed "+config_name+".csv"
    if diagnostic_level >= 3:
        print "Parsed all Config Files."

def parse_data_info(data_config_file,data_line_char,diagnostic_level):
    data_info = read_data_config(data_config_file+'.csv',data_line_char,diagnostic_level)
    pickle_object(data_config_file,data_info,diagnostic_level)
    if diagnostic_level >= 3:
        print "Parsed "+data_config_file+".csv"

def make_scatter_dict(diagnostic_level):
    quantities = read_pickle('quantities_object.pkl',diagnostic_level)
    scatter_data_dict = {}
    for quantity_key in quantities.keys():
        scatter_data_dict[int(quantity_key)] = {}
    pickle_object('scatter_data_dict',scatter_data_dict,diagnostic_level)
    if diagnostic_level >= 3:
        print"Built and Pickled Scatter Data Dictionary Object."

def add_group_data_to_scatter_data_dict(quantity_ID,group_ID,data_set,diagnostic_level):
    scatter_data_dict = read_pickle('scatter_data_dict_object.pkl',diagnostic_level)
    if int(group_ID) in scatter_data_dict[int(quantity_ID)]:
        scatter_data_dict[int(quantity_ID)][int(group_ID)].append(data_set)
    else:
        scatter_data_dict[int(quantity_ID)][int(group_ID)] = []
        scatter_data_dict[int(quantity_ID)][int(group_ID)].append(data_set)
    pickle_object('scatter_data_dict',scatter_data_dict,diagnostic_level)
    if diagnostic_level >= 3:
        print "Added Group Information to Scatter Data Dictionary Object."

def process_csv_data_file(file_name,column_names,titles_row_number,data_row_number,ind_scale,dep_scale,diagnostic_level):
    col_data = []
    fh = open(file_name, "rb")
    reader = csv.reader(fh)
    for i in range(titles_row_number-1):
        reader.next()
    titles = reader.next()
    #Strip Whitespace from ends of Column Names.
    titles = [title.strip() for title in titles]
    for i in range((data_row_number-titles_row_number)-1):
        reader.next()
    reader = csv.DictReader(fh, titles)
    data = [row for row in reader]
    
    # Mask blank values with -9999
    for row in data:
        for name in row.keys():
            if row[name] == '' or row[name] == 'NaN' or row[name] == 'inf' or row[name] == '-inf':
                row[name] = -9999.0
            else:
                row[name] = row[name]
                
    for name in column_names:
        col_data.append([float(row[name]) for row in data])
        
    fh.close()
    
    if diagnostic_level >= 3:
        print "Processed .csv File:",file_name
    
    return col_data
    

def collect_data_sets(data_info,data_directory,diagnostic_level):
    count = 0
    number_of_columns = 0
    data_sets = {}
    #import progress
    
    for key in data_info.keys():
        ind_scale = float(data_info[key]['Scale_Ind'])
        dep_scale = float(data_info[key]['Scale_Dep'])
        
        d1_column_names = []
        if data_info[key]['d1_Dep_Col_Name'][0].strip() == '[':
            d1_column_names.append(data_info[key]['d1_Ind_Col_Name'].strip())
            columns_list = eval(data_info[key]['d1_Dep_Col_Name'].strip())
            for column_name in columns_list:
                d1_column_names.append(column_name.strip())
        else:
            d1_column_names.append(data_info[key]['d1_Ind_Col_Name'].strip())
            d1_column_names.append(data_info[key]['d1_Dep_Col_Name'].strip())
        
        d1_data = process_csv_data_file(data_directory+data_info[key]['d1_Filename'],d1_column_names,int(data_info[key]['d1_Col_Name_Row']),int(data_info[key]['d1_Data_Row']),ind_scale,dep_scale,diagnostic_level)
        
        d2_column_names = []
        #print data_info[key]['d2_Filename']
        if data_info[key]['d2_Dep_Col_Name'][0].strip() == '[':
            d2_column_names.append(data_info[key]['d2_Ind_Col_Name'].strip())
            columns_list = eval(data_info[key]['d2_Dep_Col_Name'].strip())
            for column_name in columns_list:
                d2_column_names.append(column_name.strip())
        else:
            d2_column_names.append(data_info[key]['d2_Ind_Col_Name'].strip())
            d2_column_names.append(data_info[key]['d2_Dep_Col_Name'].strip())
        
        d2_data = process_csv_data_file(data_directory+data_info[key]['d2_Filename'],d2_column_names,int(data_info[key]['d2_Col_Name_Row']),int(data_info[key]['d2_Data_Row']),ind_scale,dep_scale,diagnostic_level)
        
        # Scale Data Values
        #  Scale Independent Data Values
        d1_data[0] = [value/ind_scale for value in d1_data[0]]
        d2_data[0] = [value/ind_scale for value in d2_data[0]]
        
        #  Scale Dependent Data Values
        column_count = 1
        if len(d1_column_names) <= 2:
            d1_data[1] = [value/dep_scale for value in d1_data[1]]
        else:
            for data_column in d1_data[1:]:
                d1_data[column_count] = [value/dep_scale for value in data_column]
                column_count += 1
                
        column_count = 1
        if len(d2_column_names) <= 2:
            d2_data[1] = [value/dep_scale for value in d2_data[1]]
        else:
            for data_column in d2_data[1:]:
                d2_data[column_count] = [value/dep_scale for value in data_column]
                column_count += 1
        
        data_sets[key] = [d1_data,d2_data]
        number_of_columns += ((len(d1_data)-1)+(len(d2_data)-1))
        count += 1
        
    if diagnostic_level >= 1:
        print count,"Data Info records processed."
        print "Found",number_of_columns,"columns of Y-Axis data."
    
    return data_sets

def find_start_stop_index(data_set,start_value,stop_value,start_comp,stop_comp,diagnostic_level):
    #This function is used to find index numbers for start and stop points in plotting and metric values.
    found_start = False
    found_stop = False
    found_comp_start = False
    found_comp_stop = False
    
    if diagnostic_level >= 2:
        print "\n*** Find Start/End Indexes for Independent Data Column ***"
    if diagnostic_level >= 3:
        print "Start Value:", start_value
        print "Stop Value:", stop_value
        print "Metric Start Value:", start_value
        print "Metric Stop Value:", stop_value
        
    ## Find Start Value Index Number    
    rowcounter = 0
    for value in data_set:
        if found_start == False or found_stop == False or found_comp_start == False or found_comp_stop == False:
            if value >= float(start_value) and found_start == False:
                if diagnostic_level >= 3:
                    print "Independent Data Column Starts at row #:",rowcounter
                    print "With a value of:",value
                ind_col_start_index = rowcounter
                found_start = True
            if value >= float(start_comp) and found_comp_start == False:
                if diagnostic_level >= 3:
                    print "Metric Data starts at Index #:", rowcounter, "with a value of:", value
                metric_start_index = rowcounter
                found_comp_start = True
            if float(data_set[len(data_set)-1]) <= float(stop_comp)and found_comp_stop == False:
                if diagnostic_level >= 3:
                    print "Specified end of Metric data is greater than or equal to the last value in the Independent Data Column.\nUsing last value in the Independent Data Column."
                    print "Metric Data stops at Index #:", len(data_set)-1, "with a value of:", data_set[len(data_set)-1]
                metric_end_index = len(data_set)-1
                found_comp_stop = True
            if value > float(stop_comp) and found_comp_stop == False:
                if diagnostic_level >= 3:
                    print "Metric Data stops at Index #:", rowcounter-1, "with a value of: ", data_set[rowcounter-1]
                metric_end_index = rowcounter-1
                found_comp_stop = True
            if float(data_set[len(data_set)-1]) <= float(stop_value) and found_stop == False:
                if diagnostic_level >= 3:
                    print "Specified end of Data is greater than or equal to the last value in the Independent Data Column. \nUsing last value in the Independent Data Column."
                    print "Value used is: "+str(data_set[len(data_set)-1])+"\n"
                ind_col_end_index = len(data_set)-1
                found_stop = True
            if value > float(stop_value) and found_stop == False:
                if diagnostic_level >= 3:
                    print "Independent Data Column Ends at Index #:", rowcounter, "with a value of:", data_set[rowcounter]
                ind_col_end_index = rowcounter
                found_stop = True
            rowcounter += 1
        else:
            if diagnostic_level >= 2:
                print "*** Start/End Indexes Found ***"
                print "Independent Data Col Start Index:", ind_col_start_index
                print "Independent Data Col End Index:", ind_col_end_index
                print "Metric Start Index:", metric_start_index
                print "Metric End Index:", metric_end_index
            break
    return (ind_col_start_index, ind_col_end_index, metric_start_index, metric_end_index)

def assemble_data_to_plot(diagnostic_level):
    pass


# Test Functions
def test_add_group_data_to_scatter_dict():
    add_group_data_to_scatter_data_dict(1,4,[['x1','y1'],['x2','y2']])
    add_group_data_to_scatter_data_dict(1,4,[['x3','y3'],['x4','y4']])
    
    scatter_data_dict = read_pickle('scatter_data_dict_object.pkl')
    for quantity in scatter_data_dict.keys():
        for group in scatter_data_dict[quantity].keys():
            for data_pair in scatter_data_dict[quantity][group]:
                data_x = data_pair[0]
                data_y = data_pair[1]
                print data_x,",",data_y

def test_pickle():
    style_data = read_pickle('styles_object.pkl')
    print "There are",len(style_data),"styles."
    
    group_data = read_pickle('groups_object.pkl')
    print "There are",len(group_data),"groups."
    
    quantity_data = read_pickle('quantities_object.pkl')
    print "There are",len(quantity_data),"quantities."
    
    data_info = read_pickle(data_config_file+'_object.pkl')
    print "There are",len(data_info),"data records to process."
    
    #print data_info['1~4~FM_SNL_05~T_upper~Time']['Quantity']
    #print quantity_data['1']
    #print style_data[group_data['1']['Style_ID']]['Edge_Color']
    
    #for key in sorted(data_info.keys()):
    #    print key
    

def test_process_csv_data_file():
    print "Testing the CSV file import function."
    
    filename = "/Users/bwklein/Development/FDS-SMV/Validation/FM_SNL/FDS_Output_Files/FM_SNL_21_v5_devc.csv"
    #filename = "/Users/bwklein/Development/FDS-SMV/Validation/Steckler_Compartment/Experimental_Data/Steckler_Test_10.csv"
    column_names = ['FDS Time','Sector3 Ch15','Sector2 Ch6']
    #column_names = ['Height','V_C','T_in','V_R1']
    titles_row_number = 2
    
    col_data_set = process_csv_data_file(filename,column_names,titles_row_number)
    #First list is X Axis Data
    x_data = col_data_set.pop(0)
    print "X Axis Data:",x_data
    
    #Subsequent Lists are Y Axis Data sets that pair with the X_Axis in the first list.
    print "Y Axis Data:"
    for data_col in col_data_set:
        print data_col
        print "***"


# Call Functions for Testing.
#parse_configs(config_files)
#parse_data_info(data_config_file)
#make_scatter_dict()
#test_add_group_data_to_scatter_dict()
#test_process_csv_data_file()
#test_pickle()