#!/bin/csh -f
setenv SVNROOT ~/FDS-SMV
set DDIR=$SVNROOT/Validation/Backward_Facing_Step/Current_Results
set WDIR=$SVNROOT/Validation/Backward_Facing_Step/FDS_Output_Files
cp $DDIR/*devc.csv $WDIR
cp $DDIR/*line.csv $WDIR
cp $DDIR/*svn.txt $WDIR

