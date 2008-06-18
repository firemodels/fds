import popen2
import csv
import os

config_file_name = "fds2ascii_Config_File.csv"
#path_fds2ascii = "/Applications/NIST/FDS/Utilities/fds2ascii"

#Open Config File Object
try:
    fh = file(config_file_name, 'U')
except:
    print"!!! The Config File "+config_file_name+" does not exist or the path defined in the script is incorrect. !!!"
    exit()

#Read file with csv module.
data_array = csv.reader(fh)
#Convert into List object
config_lists = [list(sublist) for sublist in data_array]
print str(len(config_lists))+" lines read in from "+config_file_name+"\n"

print config_lists

def extract_ascii_data(input_data):
    # change to working directory, the first string in each row of the config file.
    print input_data[0]
    os.chdir(input_data[0])
    print "Working Directory:",os.getcwd()
    # Open fds2ascii and create a seperate output and input stream.
    [theoutput, theinput] = popen2.popen2('fds2ascii')
    # Loop through input fields (columns) for the line passed in from the fds2ascii_Config_File.csv file.
    for field in input_data[1:]:
        if field != '':
            theinput.write(field+"\n")
            #print theoutput.readlines()
            print "Input:",field
        else:
            pass

# Pass each line in the config file to the extraction function.
for data_row in config_lists:
    extract_ascii_data(data_row)

print "*** Finished ***"