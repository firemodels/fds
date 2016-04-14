#!/bin/csh -f
 
set WDIR=$FDSSMV/Validation/Monin_Obukhov/FDS_Output_Files
set DDIR=$FDSSMV/Validation/Monin_Obukhov/Current_Results
cd $DDIR
cp *line.csv $WDIR
cp *git.txt $WDIR

