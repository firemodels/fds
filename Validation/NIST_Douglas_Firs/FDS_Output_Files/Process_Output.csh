#!/bin/csh -f
setenv SVNROOT ~/FDS-SMV
set DDIR=$SVNROOT/Validation/NIST_Douglas_Firs/Current_Results
set WDIR=$SVNROOT/Validation/NIST_Douglas_Firs/FDS_Output_Files
cp $DDIR/*hrr.csv $WDIR
cp $DDIR/*svn.txt $WDIR

