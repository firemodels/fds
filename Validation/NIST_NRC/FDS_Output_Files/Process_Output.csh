#!/bin/csh -f
setenv SVNROOT ~/FDS-SMV
set PDIR=$SVNROOT/Utilities/Data_Processing
set WDIR=$SVNROOT/Validation/NIST_NRC/FDS_Output_Files
set DDIR=$SVNROOT/Validation/NIST_NRC/Current_Results
cd $DDIR
$PDIR/layer_height < $WDIR/NIST_NRC_01_v5_HGL.input
$PDIR/layer_height < $WDIR/NIST_NRC_02_v5_HGL.input
$PDIR/layer_height < $WDIR/NIST_NRC_03_v5_HGL.input
$PDIR/layer_height < $WDIR/NIST_NRC_04_v5_HGL.input
$PDIR/layer_height < $WDIR/NIST_NRC_05_v5_HGL.input
$PDIR/layer_height < $WDIR/NIST_NRC_07_v5_HGL.input
$PDIR/layer_height < $WDIR/NIST_NRC_08_v5_HGL.input
$PDIR/layer_height < $WDIR/NIST_NRC_09_v5_HGL.input
$PDIR/layer_height < $WDIR/NIST_NRC_10_v5_HGL.input
$PDIR/layer_height < $WDIR/NIST_NRC_13_v5_HGL.input
$PDIR/layer_height < $WDIR/NIST_NRC_14_v5_HGL.input
$PDIR/layer_height < $WDIR/NIST_NRC_15_v5_HGL.input
$PDIR/layer_height < $WDIR/NIST_NRC_16_v5_HGL.input
$PDIR/layer_height < $WDIR/NIST_NRC_17_v5_HGL.input
$PDIR/layer_height < $WDIR/NIST_NRC_18_v5_HGL.input
cp NIST*HGL.csv $WDIR
cp NIST*devc.csv $WDIR
cd $WDIR


