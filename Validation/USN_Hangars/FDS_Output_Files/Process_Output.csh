#!/bin/csh -f
setenv SVNROOT ~/FDS-SMV
set DDIR=$SVNROOT/Validation/USN_Hangars/Current_Results
set WDIR=$SVNROOT/Validation/USN_Hangars/FDS_Output_Files
cp $DDIR/*devc.csv $WDIR
cp $DDIR/*svn.txt $WDIR

