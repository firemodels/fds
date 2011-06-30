#!/bin/csh -f
setenv SVNROOT ~/FDS-SMV
set PDIR=$SVNROOT/Utilities/Data_Processing
set WDIR=$SVNROOT/Validation/Hamins_CH4/FDS_Output_Files
set DDIR=$SVNROOT/Validation/Hamins_CH4/Current_Results
cd $DDIR
$PDIR/Hamins_CH4_Average
cp Hamins*devc_avg.csv $WDIR
cp *svn.txt $WDIR
cd $WDIR

