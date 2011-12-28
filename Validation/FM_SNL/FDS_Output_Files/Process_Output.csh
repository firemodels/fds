#!/bin/csh -f
setenv SVNROOT ~/FDS-SMV
set PDIR=$SVNROOT/Utilities/Data_Processing
set WDIR=$SVNROOT/Validation/FM_SNL/FDS_Output_Files
set DDIR=$SVNROOT/Validation/FM_SNL/Current_Results
cd $DDIR
$PDIR/layer_height < $WDIR/FM_SNL_01.input
$PDIR/layer_height < $WDIR/FM_SNL_02.input
$PDIR/layer_height < $WDIR/FM_SNL_03.input
$PDIR/layer_height < $WDIR/FM_SNL_04.input
$PDIR/layer_height < $WDIR/FM_SNL_05.input
$PDIR/layer_height < $WDIR/FM_SNL_06.input
$PDIR/layer_height < $WDIR/FM_SNL_07.input
$PDIR/layer_height < $WDIR/FM_SNL_08.input
$PDIR/layer_height < $WDIR/FM_SNL_09.input
$PDIR/layer_height < $WDIR/FM_SNL_10.input
$PDIR/layer_height < $WDIR/FM_SNL_11.input
$PDIR/layer_height < $WDIR/FM_SNL_12.input
$PDIR/layer_height < $WDIR/FM_SNL_13.input
$PDIR/layer_height < $WDIR/FM_SNL_14.input
$PDIR/layer_height < $WDIR/FM_SNL_15.input
$PDIR/layer_height < $WDIR/FM_SNL_16.input
$PDIR/layer_height < $WDIR/FM_SNL_17.input
$PDIR/layer_height < $WDIR/FM_SNL_21.input
$PDIR/layer_height < $WDIR/FM_SNL_22.input
cp FM_SNL*HGL.csv $WDIR
cp FM_SNL*devc.csv $WDIR
cp *svn.txt $WDIR
cd $WDIR


