#!/bin/csh -f
setenv SVNROOT ~/FDS-SMV
set DDIR=$SVNROOT/Validation/Purdue_Flames/Current_Results
set WDIR=$SVNROOT/Validation/Purdue_Flames/FDS_Output_Files
#cp $DDIR/*devc.csv $WDIR
cp $DDIR/*line.csv $WDIR
cp $DDIR/*svn.txt $WDIR

