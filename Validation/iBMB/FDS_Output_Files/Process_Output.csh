#!/bin/csh -f
setenv SVNROOT ~/FDS-SMV
set PDIR=$SVNROOT/Utilities/Data_Processing
set WDIR=$SVNROOT/Validation/iBMB/FDS_Output_Files
set DDIR=$SVNROOT/Validation/iBMB/Current_Results
cd $DDIR
$PDIR/layer_height < $WDIR/ICFMP4_01_HGL.input
$PDIR/layer_height < $WDIR/ICFMP5_04_HGL.input
cp *HGL.csv  $WDIR
cp *devc.csv $WDIR
cp *svn.txt $WDIR
cd $WDIR

