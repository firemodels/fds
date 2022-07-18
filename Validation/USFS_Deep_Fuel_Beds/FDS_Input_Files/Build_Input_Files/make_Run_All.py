# 7/18/2022 Noelle Crump, to generate Run_All.sh 
# this one is for USFS Deep Fuel Beds

import shutil
import os
from os import listdir

input_file_list = [f for f in listdir('../') if f[-4:] == '.fds']
shutil.copy('Run_All_base.txt','Run_All.sh')
forestring = '$QFDS $DEBUG $QUEUE -p 16 -d $INDIR '

f = open('Run_All.sh', 'a')
for filename in input_file_list:
    f.write(forestring + filename +'\n')
f.write('\necho FDS cases submitted')
f.close()

# Move Run_All.sh  files up to Case Folder
os.replace("Run_All.sh", "../../Run_All.sh")
os.system('chmod +x  ../../Run_All.sh')