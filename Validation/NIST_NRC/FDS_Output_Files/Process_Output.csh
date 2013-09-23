#!/bin/csh -f
setenv SVNROOT ~/FDS-SMV
set DDIR=$SVNROOT/Validation/NIST_NRC/Current_Results
set WDIR=$SVNROOT/Validation/NIST_NRC/FDS_Output_Files
cp $DDIR/*devc.csv $WDIR
cp $DDIR/*svn.txt $WDIR

