#!/bin/csh -f
set PDIR=~/FDS_Repository/Utilities/Data_Processing
set WDIR=~/FDS_Repository/Validation/VTT/FDS_Output_Files
set DDIR=~/VALIDATION/VTT/FDS_5.3
cd $DDIR
$PDIR/layer_height < $WDIR/VTT_01_v5_HGL.input
$PDIR/layer_height < $WDIR/VTT_02_v5_HGL.input
$PDIR/layer_height < $WDIR/VTT_03_v5_HGL.input
cp VTT*HGL.csv $WDIR
cp VTT*devc.csv $WDIR
cd $WDIR


