#!/bin/csh -f
set PDIR=~/FDS_Repository/Utilities/Data_Processing
set WDIR=~/FDS_Repository/Validation/NIST_NRC/FDS_Output_Files
set DDIR=~/VALIDATION/NIST_NRC/FDS_5.3
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


