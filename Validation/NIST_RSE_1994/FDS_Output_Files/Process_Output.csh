#!/bin/csh -f
setenv SVNROOT ~/FDS-SMV
set WDIR=$SVNROOT/Validation/NIST_RSE_1994/FDS_Output_Files
set DDIR=$SVNROOT/Validation/NIST_RSE_1994/Current_Results
cp $DDIR/*devc.csv $WDIR
cp $DDIR/*svn.txt $WDIR

