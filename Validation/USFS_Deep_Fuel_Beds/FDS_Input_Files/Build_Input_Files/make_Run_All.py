# 7/18/2022 Noelle Crump, to generate Run_All.sh 
# this one is for USFS Deep Fuel Beds

import shutil
import os
from os import listdir

input_file_list = [f for f in listdir('../') if f[-4:] == '.fds']
forestring = '$QFDS $DEBUG $QUEUE -p 16 -d $INDIR '

f = open('Run_All.sh', 'w')
f.write('#!/bin/bash\n\n')
f.write('# This script runs a set of Validation Cases on a Linux machine with a batch queuing system.\n')
f.write('# See the file Validation/Common_Run_All.sh for more information.\n')
f.write('export SVNROOT=`pwd`/../..\nsource $SVNROOT/Validation/Common_Run_All.sh\n\n')
f.write('# good burn: 136, 114, marginal burn: 41,80 no burn: 30,60)\n\n')

for filename in input_file_list:
    f.write(forestring + filename +'\n')
f.write('\necho FDS cases submitted')
f.close()

# Move Run_All.sh  files up to Case Folder
os.replace("Run_All.sh", "../../Run_All.sh")
os.system('chmod +x  ../../Run_All.sh')