#!/bin/csh -f
setenv SVNROOT ~/FDS-SMV
set DDIR=$SVNROOT/Validation/ATF_Corridors/Current_Results
set WDIR=$SVNROOT/Validation/ATF_Corridors/FDS_Output_Files
cp $DDIR/ATF*devc.csv $WDIR

