#!/bin/csh -f
 
set DDIR=$FDSSMV/Validation/DelCo_Trainers/Current_Results
set WDIR=$FDSSMV/Validation/DelCo_Trainers/FDS_Output_Files
cp $DDIR/*devc.csv $WDIR
cp $DDIR/*git.txt $WDIR

