#!/bin/csh -f
setenv SVNROOT ~/FDS-SMV
set WDIR=$SVNROOT/Validation/Heskestad_Flame_Height/FDS_Output_Files
set DDIR=$SVNROOT/Validation/Heskestad_Flame_Height/Current_Results
cd $DDIR
cp Qs*line.csv $WDIR
cp Qs*hrr.csv  $WDIR
cp *svn.txt    $WDIR
cd $WDIR
cd ./

