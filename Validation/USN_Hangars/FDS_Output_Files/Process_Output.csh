#!/bin/csh -f
 
set DDIR=$FDSSMV/Validation/USN_Hangars/Current_Results
set WDIR=$FDSSMV/Validation/USN_Hangars/FDS_Output_Files
cp $DDIR/*devc.csv $WDIR
cp $DDIR/*git.txt $WDIR

