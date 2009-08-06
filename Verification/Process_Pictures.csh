#!/bin/csh -f
setenv SVNROOT ~/FDS-SMV
set RUNSMV=$SVNROOT/SMV_5/Build/INTEL_LINUX_32/smokeview
cd Fires
$RUNSMV -runscript room_fire
cd ..
cd NS_Analytical_Solution
$RUNSMV -runscript ns2d_64
cd ..


