#!/bin/csh -f
setenv SVNROOT ~/FDS-SMV
set PDIR=$SVNROOT/Utilities/Data_Processing
set WDIR=$SVNROOT/Validation/FM_SNL/FDS_Output_Files
set DDIR=$SVNROOT/Validation/FM_SNL/Current_Results
cd $DDIR
$PDIR/layer_height < $WDIR/FM_SNL_04.input
$PDIR/layer_height < $WDIR/FM_SNL_05.input
$PDIR/layer_height < $WDIR/FM_SNL_21.input
cp FM_SNL*HGL.csv $WDIR
cp FM_SNL*devc.csv $WDIR
cd $WDIR


