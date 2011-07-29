#!/bin/csh -f
setenv SVNROOT ~/FDS-SMV
set PDIR=$SVNROOT/Utilities/Data_Processing
set WDIR=$SVNROOT/Validation/Harrison_Spill_Plumes/FDS_Output_Files
set DDIR=$SVNROOT/Validation/Harrison_Spill_Plumes/Current_Results
cd $DDIR
cp *devc.csv  $WDIR
cp *svn.txt   $WDIR
cd $WDIR


