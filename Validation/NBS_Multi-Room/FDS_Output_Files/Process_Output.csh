#!/bin/csh -f
setenv SVNROOT ~/FDS-SMV
set DDIR=$SVNROOT/Validation/NBS_Multi-Room/Current_Results
set WDIR=$SVNROOT/Validation/NBS_Multi-Room/FDS_Output_Files
cp $DDIR/NBS*devc.csv $WDIR
cp $DDIR/*svn.txt $WDIR

