#!/bin/csh -f
 
set WDIR=$FDSSMV/Validation/Restivo_Experiment/FDS_Output_Files
set DDIR=$FDSSMV/Validation/Restivo_Experiment/Current_Results
cd $DDIR
cp Res*line.csv $WDIR
cp *git.txt $WDIR

