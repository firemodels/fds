#!/bin/csh -f
setenv SVNROOT ~/FDS-SMV
set PDIR=$SVNROOT/Utilities/Data_Processing
set WDIR=$SVNROOT/Validation/ATF_Corridors/FDS_Output_Files
set DDIR=$SVNROOT/Validation/ATF_Corridors/Current_Results
cd $DDIR
$PDIR/layer_height < $WDIR/ATF_Corridors_050_kW_HGL_1.input
$PDIR/layer_height < $WDIR/ATF_Corridors_050_kW_HGL_2.input
$PDIR/layer_height < $WDIR/ATF_Corridors_100_kW_HGL_1.input
$PDIR/layer_height < $WDIR/ATF_Corridors_100_kW_HGL_2.input
$PDIR/layer_height < $WDIR/ATF_Corridors_240_kW_HGL_1.input
$PDIR/layer_height < $WDIR/ATF_Corridors_240_kW_HGL_2.input
$PDIR/layer_height < $WDIR/ATF_Corridors_250_kW_HGL_1.input
$PDIR/layer_height < $WDIR/ATF_Corridors_250_kW_HGL_2.input
$PDIR/layer_height < $WDIR/ATF_Corridors_500_kW_HGL_1.input
$PDIR/layer_height < $WDIR/ATF_Corridors_500_kW_HGL_2.input
$PDIR/layer_height < $WDIR/ATF_Corridors_Mix_kW_HGL_1.input
$PDIR/layer_height < $WDIR/ATF_Corridors_Mix_kW_HGL_2.input
cp ATF*HGL.csv $WDIR
cp ATF*devc.csv $WDIR
cp *svn.txt $WDIR
cd $WDIR


