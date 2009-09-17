#!/bin/csh -f
setenv SVNROOT ~/FDS-SMV
setenv FDS $SVNROOT/FDS_Compilation/intel_linux_32/fds5_intel_linux_32
setenv RUNFDS $SVNROOT/Utilities/Scripts/runfds.csh
setenv BASEDIR `pwd`

./FDS_Cases.csh

echo FDS cases submitted
