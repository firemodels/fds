#!/bin/csh -f

#  This script runs all FDS and Smokeview verifciation cases.
#  If you decide you want to stop all the runs before they would
#  finish normally then type:
#
#  setenv STOPFDS
#
#  at a command prompt then re-run this script. After FDS jobs have stopped
#  then open up a NEW shell and re-run script.
#
setenv SVNROOT ~/FDS-SMV
setenv FDS $SVNROOT/FDS_Compilation/intel_linux_32/fds5_intel_linux_32
setenv RUNFDS $SVNROOT/Utilities/Scripts/runfds.csh
setenv BASEDIR `pwd`

echo FDS cases submitted

./FDS_Cases.csh

cd $BASEDIR/SMV_Scripts

./SMV_Cases.csh
