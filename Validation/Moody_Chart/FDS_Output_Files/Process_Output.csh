#!/bin/csh -f
 
set DDIR=$FDSSMV/Validation/Moody_Chart/Current_Results
set WDIR=$FDSSMV/Validation/Moody_Chart/FDS_Output_Files
cp $DDIR/*devc.csv $WDIR
cp $DDIR/*git.txt $WDIR

