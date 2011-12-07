#!/bin/csh -f
setenv SVNROOT ~/FDS-SMV
set WDIR=$SVNROOT/Validation/Lattimer_Tilted_Wall/FDS_Output_Files
set DDIR=$SVNROOT/Validation/Lattimer_Tilted_Wall/Current_Results
cd $DDIR
cp *devc.csv $WDIR
cp *svn.txt $WDIR

