#!/bin/csh -f
 
set WDIR=$FDSSMV/Validation/Bryant_Doorway/FDS_Output_Files
set DDIR=$FDSSMV/Validation/Bryant_Doorway/Current_Results
cp $DDIR/Br*line.csv $WDIR
cp $DDIR/*git.txt $WDIR


