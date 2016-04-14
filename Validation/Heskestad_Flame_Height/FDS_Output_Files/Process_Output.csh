#!/bin/csh -f
 
set WDIR=$FDSSMV/Validation/Heskestad_Flame_Height/FDS_Output_Files
set DDIR=$FDSSMV/Validation/Heskestad_Flame_Height/Current_Results
cd $DDIR
cp Qs*line.csv $WDIR
cp Qs*hrr.csv  $WDIR
cp *git.txt    $WDIR
cd $WDIR
cd ./

