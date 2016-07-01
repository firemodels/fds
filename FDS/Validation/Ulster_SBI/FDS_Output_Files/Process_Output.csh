#!/bin/csh -f
 
set WDIR=$FDSSMV/Validation/Ulster_SBI/FDS_Output_Files
set DDIR=$FDSSMV/Validation/Ulster_SBI/Current_Results
cd $DDIR
cp Ulster*line.csv $WDIR
cp *git.txt $WDIR

