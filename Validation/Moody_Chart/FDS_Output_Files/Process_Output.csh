#!/bin/csh -f
setenv SVNROOT ~/FDS-SMV
set DDIR=$SVNROOT/Validation/Moody_Chart/Current_Results
set WDIR=$SVNROOT/Validation/Moody_Chart/FDS_Output_Files
cp $DDIR/*devc.csv $WDIR

