#!/bin/csh -f
 
set WDIR=$FDSSMV/Validation/FM_SNL/FDS_Output_Files
set DDIR=$FDSSMV/Validation/FM_SNL/Current_Results
cp $DDIR/*devc.csv $WDIR
cp $DDIR/*git.txt $WDIR


