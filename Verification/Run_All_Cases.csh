#!/bin/csh -f

#  This script runs all FDS and Smokeview verifciation cases.
#
#  To abort FDS runs, uncomment the following line and rerun this script.
#
#setenv STOPFDS
#
setenv SVNROOT ~/FDS-SMV
setenv FDS $SVNROOT/FDS_Compilation/intel_linux_32/fds5_intel_linux_32
setenv SMOKEZIP $SVNROOT/Utilities/smokezip/INTEL_LINUX_32/smokezip_linux_32
setenv BASEDIR `pwd`

#setenv RUNFDS $SVNROOT/Utilities/Scripts/runfds.csh

# To run the "background" version of this script, 
# comment the line above contaiing RUNFDS
# uncomment the two lines below containing RUNFDS and BACKGROUND

setenv RUNFDS $SVNROOT/Utilities/Scripts/runfds_bg.csh
setenv BACKGROUND "$SVNROOT/Utilities/background/INTEL_LINUX_32/background -hosts ~/run_nodes -u 75 -d 5"

setenv RUNFDSFG $SVNROOT/Utilities/Scripts/runfds_fg.csh

echo FDS cases submitted

./FDS_Cases.csh

cd $BASEDIR/scripts

./SMV_Cases.csh

cd $BASEDIR/scripts
./run_wui_tree_test.csh
