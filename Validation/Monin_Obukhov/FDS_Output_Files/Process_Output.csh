#!/bin/csh -f
setenv SVNROOT ~/FDS-SMV
set WDIR=$SVNROOT/Validation/Monin_Obukhov/FDS_Output_Files
set DDIR=$SVNROOT/Validation/Monin_Obukhov/Current_Results
cd $DDIR
cp *line.csv $WDIR
cp *svn.txt $WDIR

