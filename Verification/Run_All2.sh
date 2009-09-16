#!/bin/bash
export SVNROOT=~/FDS-SMV
export FDS=$SVNROOT/FDS_Compilation/intel_linux_64/fds5_intel_linux_64
export RUNFDS=$SVNROOT/Utilities/Scripts/runfds.sh
export BASEDIR=`pwd`

DELIM='/'

FdsFileDirs=`ls -1 -d */ | cut -d'/' -f 1`

for dir_now in ${FdsFileDirs}
   do
   echo 'Directory is now' ${dir_now}
   FdsFiles=`ls -1 -d  ${dir_now}${DELIM}*.fds | cut -d'.' -f 1 | cut -d${DELIM} -f 2` 
   for fds_now in ${FdsFiles}
     do
     echo 'Fds file is now ' ${fds_now}'.fds'
     $RUNFDS ${dir_now} ${fds_now}
     done
done

echo FDS cases submitted
