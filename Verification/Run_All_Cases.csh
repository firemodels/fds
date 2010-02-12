#!/bin/csh -f
setenv SVNROOT ~/FDS-SMV
setenv FDS $SVNROOT/FDS_Compilation/intel_osx_32/fds5_intel_osx_32
setenv RUNFDS $SVNROOT/Utilities/Scripts/runfdsosx.csh
setenv BASEDIR `pwd`

echo FDS cases submitted

./FDS_Cases.csh

cd $BASEDIR

#./SMV_Cases.csh

