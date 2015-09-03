#!/bin/csh -f
 
set WDIR=$FDSSMV/Validation/NRCC_Facade/FDS_Output_Files
set DDIR=$FDSSMV/Validation/NRCC_Facade/Current_Results
cp $DDIR/NRCC*line.csv $WDIR
cp $DDIR/*git.txt $WDIR

