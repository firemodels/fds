#!/bin/csh -f
set PDIR=~/FDS_Repository/Utilities/Data_Processing
set WDIR=~/FDS_Repository/Validation/Hamins_CH4/FDS_Output_Files
set DDIR=~/FDS_Repository/Validation/Hamins_CH4/Current_Results
cd $DDIR
$PDIR/Hamins_CH4_Average
cp Hamins*devc_avg.csv $WDIR
cd $WDIR

