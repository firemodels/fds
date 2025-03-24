
import glob
import re

print("Series,Case,Wall Clock Time (h),Processes,CPU hours,Git Hash,Date")
files = glob.glob('*/Current_Results/*.out')
for file in files:
    with open(file, 'r') as f:
        for line in f:
            if re.search("Number of MPI Processes", line):
                np = line[27:30]
            if re.search("Revision         :", line):
                githash = line[20:60]
            if re.search("Revision Date", line):
                revdate = line[20:60]
            if re.search("Total Elapsed", line):
                file_name=re.sub("/Current_Results/",",",f.name)
                file_name2=re.sub(".out","",file_name)
                print(file_name2,",",float(line[36:48])/3600,",",int(np),",",int(np)*float(line[36:48])/3600,",",githash.rstrip(),",",revdate.rstrip())

