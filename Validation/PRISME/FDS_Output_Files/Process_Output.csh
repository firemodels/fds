#!/bin/csh -f
setenv SVNROOT ~/FDS-SMV
set DDIR=$SVNROOT/Validation/PRISME/Current_Results
set WDIR=$SVNROOT/Validation/PRISME/FDS_Output_Files
cp $DDIR/*devc.csv $WDIR

