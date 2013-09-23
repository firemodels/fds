#!/bin/csh -f
setenv SVNROOT ~/FDS-SMV
set PDIR=$SVNROOT/Utilities/Data_Processing
set WDIR=$SVNROOT/Validation/FM_SNL/FDS_Output_Files
set DDIR=$SVNROOT/Validation/FM_SNL/Current_Results
cd $DDIR
$PDIR/layer_height < $WDIR/FM_SNL_01_HGL.input
$PDIR/layer_height < $WDIR/FM_SNL_02_HGL.input
$PDIR/layer_height < $WDIR/FM_SNL_03_HGL.input
$PDIR/layer_height < $WDIR/FM_SNL_04_HGL.input
$PDIR/layer_height < $WDIR/FM_SNL_05_HGL.input
$PDIR/layer_height < $WDIR/FM_SNL_06_HGL.input
$PDIR/layer_height < $WDIR/FM_SNL_07_HGL.input
$PDIR/layer_height < $WDIR/FM_SNL_08_HGL.input
$PDIR/layer_height < $WDIR/FM_SNL_09_HGL.input
$PDIR/layer_height < $WDIR/FM_SNL_10_HGL.input
$PDIR/layer_height < $WDIR/FM_SNL_11_HGL.input
$PDIR/layer_height < $WDIR/FM_SNL_12_HGL.input
$PDIR/layer_height < $WDIR/FM_SNL_13_HGL.input
$PDIR/layer_height < $WDIR/FM_SNL_14_HGL.input
$PDIR/layer_height < $WDIR/FM_SNL_15_HGL.input
$PDIR/layer_height < $WDIR/FM_SNL_16_HGL.input
$PDIR/layer_height < $WDIR/FM_SNL_17_HGL.input
$PDIR/layer_height < $WDIR/FM_SNL_21_HGL.input
$PDIR/layer_height < $WDIR/FM_SNL_22_HGL.input
cp FM_SNL*HGL.csv $WDIR
cp FM_SNL*devc.csv $WDIR
cp *svn.txt $WDIR
cd $WDIR


