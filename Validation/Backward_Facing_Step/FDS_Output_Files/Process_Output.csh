#!/bin/csh -f
 
set DDIR=$FDSSMV/Validation/Backward_Facing_Step/Current_Results
set WDIR=$FDSSMV/Validation/Backward_Facing_Step/FDS_Output_Files
#cp $DDIR/*devc.csv $WDIR
cp $DDIR/*line.csv $WDIR
cp $DDIR/*git.txt $WDIR

