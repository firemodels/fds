#!/bin/csh -f
setenv SVNROOT ~/FDS-SMV
set WDIR=$SVNROOT/Validation/LLNL_Enclosure/FDS_Output_Files
set DDIR=$SVNROOT/Validation/LLNL_Enclosure/Current_Results
cd $DDIR
cp *devc.csv $WDIR
cp *svn.txt $WDIR
cd $WDIR


