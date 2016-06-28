#!/bin/csh -f
 
set WDIR=$FDSSMV/Validation/VTT/FDS_Output_Files
set DDIR=$FDSSMV/Validation/VTT/Current_Results
cp $DDIR/*devc.csv $WDIR
cp $DDIR/*git.txt $WDIR


