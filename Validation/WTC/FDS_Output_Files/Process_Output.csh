#!/bin/csh -f
setenv SVNROOT ~/FDS-SMV
set PDIR=$SVNROOT/Utilities/Data_Processing
set WDIR=$SVNROOT/Validation/WTC/FDS_Output_Files
set DDIR=$SVNROOT/Validation/WTC/Current_Results
cd $DDIR
$PDIR/layer_height < $WDIR/WTC_01_HGL.input
$PDIR/layer_height < $WDIR/WTC_02_HGL.input
$PDIR/layer_height < $WDIR/WTC_03_HGL.input
$PDIR/layer_height < $WDIR/WTC_04_HGL.input
$PDIR/layer_height < $WDIR/WTC_05_HGL.input
$PDIR/layer_height < $WDIR/WTC_06_HGL.input
cp WTC*HGL.csv $WDIR
cp WTC*devc.csv $WDIR
cd $WDIR


