#!/bin/csh -f
setenv SVNROOT ~/FDS-SMV
set WDIR=$SVNROOT/Validation/ATF_Corridors/FDS_Output_Files
set DDIR=$SVNROOT/Validation/ATF_Corridors/Current_Results
cp $DDIR/*devc.csv $WDIR
cp $DDIR/*svn.txt $WDIR


