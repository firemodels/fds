#!/bin/csh -f
 
set WDIR=$FDSSMV/Validation/Beyler_Hood/FDS_Output_Files
set DDIR=$FDSSMV/Validation/Beyler_Hood/Current_Results
cp $DDIR/*devc.csv $WDIR
cp $DDIR/*git.txt $WDIR

