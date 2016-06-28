#!/bin/csh -f
 
set DDIR=$FDSSMV/Validation/NBS_Multi-Room/Current_Results
set WDIR=$FDSSMV/Validation/NBS_Multi-Room/FDS_Output_Files
cp $DDIR/NBS*devc.csv $WDIR
cp $DDIR/*git.txt $WDIR

