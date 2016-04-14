#!/bin/csh -f
 
set DDIR=$FDSSMV/Validation/PRISME/Current_Results
set WDIR=$FDSSMV/Validation/PRISME/FDS_Output_Files
cp $DDIR/*devc.csv $WDIR
cp $DDIR/*git.txt $WDIR
cp $DDIR/PRS*hrr.csv $WDIR

