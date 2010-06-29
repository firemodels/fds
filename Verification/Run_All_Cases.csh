#!/bin/csh -f

#  This script runs all FDS and Smokeview verifciation cases.
#
#  To abort FDS runs, uncomment the following line and rerun this script.
#
#setenv STOPFDS
#
setenv SVNROOT ~/FDS-SMV
setenv FDS $SVNROOT/FDS_Compilation/intel_linux_32/fds5_intel_linux_32
setenv BASEDIR `pwd`

setenv RUNFDS $SVNROOT/Utilities/Scripts/runfds.csh

# To run the "background" version of this script, 
# uncomment the two lines below containing RUNFDS and BACKGROUND

#setenv RUNFDS $SVNROOT/Utilities/Scripts/runfds_bg.csh
#setenv BACKGROUND "$SVNROOT/Utilities/background/INTEL_LINUX_32/background -u 50 -d 10"

echo FDS cases submitted

./FDS_Cases.csh

cd $BASEDIR/SMV_Scripts

./SMV_Cases.csh
