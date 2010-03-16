#!/bin/csh -f
setenv SVNROOT ~/FDS-SMV
setenv FDS $SVNROOT/FDS_Compilation/intel_linux_64/fds5_intel_linux_64
set RUNFDS=$SVNROOT/Utilities/Scripts/runfds_sge.csh
setenv BASEDIR `pwd`

cp $BASEDIR/FDS_Input_Files/sge-fds-array.sh $BASEDIR/Current_Results

$RUNFDS Current_Results BRE_Spray_ 1-24

echo FDS cases submitted
