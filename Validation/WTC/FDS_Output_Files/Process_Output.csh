#!/bin/csh -f
set PDIR=~/FDS_Repository/Utilities/Data_Processing
set WDIR=~/FDS_Repository/Validation/WTC/FDS_Output_Files
set DDIR=~/VALIDATION/WTC/FDS_5.3
cd $DDIR
$PDIR/layer_height < $WDIR/WTC_01_v5_HGL.input
$PDIR/layer_height < $WDIR/WTC_02_v5_HGL.input
$PDIR/layer_height < $WDIR/WTC_03_v5_HGL.input
$PDIR/layer_height < $WDIR/WTC_04_v5_HGL.input
$PDIR/layer_height < $WDIR/WTC_05_v5_HGL.input
$PDIR/layer_height < $WDIR/WTC_06_v5_HGL.input
cp WTC*HGL.csv $WDIR
cp WTC*devc.csv $WDIR
cd $WDIR


