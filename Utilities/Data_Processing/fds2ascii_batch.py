import subprocess
import csv
import os

config_file_name = "fds2ascii_batch_Config_File.csv"
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

def extract_ascii_data(input_data):
    # change to working directory, the first string in each row of the config file.
    #print input_data[0]
    os.chdir(input_data[0])
    print "Working Directory:",os.getcwd()
    # Open fds2ascii and create a seperate output and input stream.
    proc = subprocess.Popen('fds2ascii',
                           shell=True,
                           stdin=subprocess.PIPE,
                           stdout=subprocess.PIPE,
                           )
    # Loop through input fields (columns) for the line passed in from the fds2ascii_Config_File.csv file.
    for field in input_data[1:]:
        if field != '':
            print "Input:",field
            proc.stdin.write('%s\n' % field)
        else:
            pass
    proc.communicate()[0]

# Pass each line in the config file to the extraction function.
for data_row in config_lists:
    extract_ascii_data(data_row)

print "*** Finished ***"
# Script Written on 06.19.08 by bryan.klein@nist.gov
# For usage instructions, see: http://code.google.com/p/fds-smv/wiki/fds2asciiBatchProcess