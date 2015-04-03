#!/bin/csh -f
setenv SVNROOT ~/FDS-SMV
set WDIR=$SVNROOT/Validation/UMD_Polymers/FDS_Output_Files
set DDIR=$SVNROOT/Validation/UMD_Polymers/Current_Results
cd $DDIR
cp *devc.csv $WDIR
cp *svn.txt $WDIR
cd $WDIR

