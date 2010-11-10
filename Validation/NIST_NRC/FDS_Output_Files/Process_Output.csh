#!/bin/csh -f
setenv SVNROOT ~/FDS-SMV
set PDIR=$SVNROOT/Utilities/Data_Processing
set WDIR=$SVNROOT/Validation/NIST_NRC/FDS_Output_Files
set DDIR=$SVNROOT/Validation/NIST_NRC/Current_Results
cd $DDIR
$PDIR/layer_height < $WDIR/NIST_NRC_01_HGL.input
$PDIR/layer_height < $WDIR/NIST_NRC_02_HGL.input
$PDIR/layer_height < $WDIR/NIST_NRC_03_HGL.input
$PDIR/layer_height < $WDIR/NIST_NRC_04_HGL.input
$PDIR/layer_height < $WDIR/NIST_NRC_05_HGL.input
$PDIR/layer_height < $WDIR/NIST_NRC_07_HGL.input
$PDIR/layer_height < $WDIR/NIST_NRC_08_HGL.input
$PDIR/layer_height < $WDIR/NIST_NRC_09_HGL.input
$PDIR/layer_height < $WDIR/NIST_NRC_10_HGL.input
$PDIR/layer_height < $WDIR/NIST_NRC_13_HGL.input
$PDIR/layer_height < $WDIR/NIST_NRC_14_HGL.input
$PDIR/layer_height < $WDIR/NIST_NRC_15_HGL.input
$PDIR/layer_height < $WDIR/NIST_NRC_16_HGL.input
$PDIR/layer_height < $WDIR/NIST_NRC_17_HGL.input
$PDIR/layer_height < $WDIR/NIST_NRC_18_HGL.input
cp NIST*HGL.csv $WDIR
cp NIST*devc.csv $WDIR
cd $WDIR


