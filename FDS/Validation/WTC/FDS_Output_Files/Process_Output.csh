#!/bin/csh -f
 
set WDIR=$FDSSMV/Validation/WTC/FDS_Output_Files
set DDIR=$FDSSMV/Validation/WTC/Current_Results
cp $DDIR/*devc.csv $WDIR
cp $DDIR/*git.txt $WDIR


