#!/bin/csh -f
setenv SVNROOT ~/FDS-SMV
set WDIR=$SVNROOT/Validation/FM_SNL/FDS_Output_Files
set DDIR=$SVNROOT/Validation/FM_SNL/Current_Results
cp $DDIR/*devc.csv $WDIR
cp $DDIR/*svn.txt $WDIR


