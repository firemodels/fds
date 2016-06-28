#!/bin/csh -f
 
set WDIR=$FDSSMV/Validation/ATF_Corridors/FDS_Output_Files
set DDIR=$FDSSMV/Validation/ATF_Corridors/Current_Results
cp $DDIR/*devc.csv $WDIR
cp $DDIR/*git.txt $WDIR


