#!/bin/csh -f
setenv SVNROOT ~/FDS-SMV
set PDIR=$SVNROOT/Utilities/Data_Processing
set WDIR=$SVNROOT/Validation/VTT/FDS_Output_Files
set DDIR=$SVNROOT/Validation/VTT/Current_Results
cd $DDIR
$PDIR/layer_height < $WDIR/VTT_01_HGL.input
$PDIR/layer_height < $WDIR/VTT_02_HGL.input
$PDIR/layer_height < $WDIR/VTT_03_HGL.input
cp VTT*HGL.csv $WDIR
cp VTT*devc.csv $WDIR
cd $WDIR


