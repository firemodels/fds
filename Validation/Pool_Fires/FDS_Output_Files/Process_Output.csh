#!/bin/csh -f
setenv SVNROOT ~/FDS-SMV
set DDIR=$SVNROOT/Validation/Pool_Fires/Current_Results
set WDIR=$SVNROOT/Validation/Pool_Fires/FDS_Output_Files
cp $DDIR/*hrr.csv $WDIR
cp $DDIR/*svn.txt $WDIR

