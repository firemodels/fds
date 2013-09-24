#!/bin/csh -f
setenv SVNROOT ~/FDS-SMV
set WDIR=$SVNROOT/Validation/UL_NIST_Vents/FDS_Output_Files
set DDIR=$SVNROOT/Validation/UL_NIST_Vents/Current_Results
cp $DDIR/*devc.csv $WDIR
cp $DDIR/*svn.txt $WDIR


