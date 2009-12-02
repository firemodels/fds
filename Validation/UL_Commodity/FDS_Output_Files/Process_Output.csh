#!/bin/csh -f
setenv SVNROOT ~/FDS-SMV
set DDIR=$SVNROOT/Validation/UL_Commodity/Current_Results
set WDIR=$SVNROOT/Validation/UL_Commodity/FDS_Output_Files
cp $DDIR/Group*hrr.csv $WDIR

