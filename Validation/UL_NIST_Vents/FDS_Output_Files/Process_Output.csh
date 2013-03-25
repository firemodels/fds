#!/bin/csh -f
setenv SVNROOT ~/FDS-SMV
set PDIR=$SVNROOT/Utilities/Data_Processing
set WDIR=$SVNROOT/Validation/UL_NIST_Vents/FDS_Output_Files
set DDIR=$SVNROOT/Validation/UL_NIST_Vents/Current_Results
cd $DDIR
$PDIR/layer_height < $WDIR/UL_NIST_Vents_Test_1_HGL.input
$PDIR/layer_height < $WDIR/UL_NIST_Vents_Test_2_HGL.input
$PDIR/layer_height < $WDIR/UL_NIST_Vents_Test_3_HGL.input
$PDIR/layer_height < $WDIR/UL_NIST_Vents_Test_4_HGL.input
cp UL*HGL.csv $WDIR
cp UL*devc.csv $WDIR
cp *svn.txt $WDIR
cd $WDIR


