#!/bin/csh -f
 
set DDIR=$FDSSMV/Validation/CAROLFIRE/Current_Results
set WDIR=$FDSSMV/Validation/CAROLFIRE/FDS_Output_Files
cp $DDIR/CAROLFIRE*devc.csv $WDIR
cp $DDIR/*git.txt $WDIR

