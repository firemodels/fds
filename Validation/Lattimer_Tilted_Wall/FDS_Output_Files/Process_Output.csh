#!/bin/csh -f
 
set WDIR=$FDSSMV/Validation/Lattimer_Tilted_Wall/FDS_Output_Files
set DDIR=$FDSSMV/Validation/Lattimer_Tilted_Wall/Current_Results
cd $DDIR
cp *line.csv $WDIR
cp *devc.csv $WDIR
cp *git.txt $WDIR

