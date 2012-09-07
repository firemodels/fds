#!/bin/csh -f
setenv SVNROOT ~/FDS-SMV
set WDIR=$SVNROOT/Validation/Beyler_Hood/FDS_Output_Files
set DDIR=$SVNROOT/Validation/Beyler_Hood/Current_Results
cp $DDIR/*devc.csv $WDIR
cp $DDIR/*svn.txt $WDIR

