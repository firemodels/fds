#!/bin/csh -f
setenv SVNROOT ~/FDS-SMV
set PDIR=$SVNROOT/Utilities/Data_Processing
set WDIR=$SVNROOT/Validation/ATF_Corridors/FDS_Output_Files
set DDIR=$SVNROOT/Validation/ATF_Corridors/Current_Results
cd $DDIR
$PDIR/layer_height < $WDIR/ATF_Corridors_050_kW_1_HGL.input
$PDIR/layer_height < $WDIR/ATF_Corridors_050_kW_2_HGL.input
$PDIR/layer_height < $WDIR/ATF_Corridors_100_kW_1_HGL.input
$PDIR/layer_height < $WDIR/ATF_Corridors_100_kW_2_HGL.input
$PDIR/layer_height < $WDIR/ATF_Corridors_240_kW_1_HGL.input
$PDIR/layer_height < $WDIR/ATF_Corridors_240_kW_2_HGL.input
$PDIR/layer_height < $WDIR/ATF_Corridors_250_kW_1_HGL.input
$PDIR/layer_height < $WDIR/ATF_Corridors_250_kW_2_HGL.input
$PDIR/layer_height < $WDIR/ATF_Corridors_500_kW_1_HGL.input
$PDIR/layer_height < $WDIR/ATF_Corridors_500_kW_2_HGL.input
$PDIR/layer_height < $WDIR/ATF_Corridors_Mix_kW_1_HGL.input
$PDIR/layer_height < $WDIR/ATF_Corridors_Mix_kW_2_HGL.input
cp ATF*HGL_1.csv $WDIR
cp ATF*HGL_2.csv $WDIR
cp ATF*devc.csv $WDIR
cp *svn.txt $WDIR
cd $WDIR


