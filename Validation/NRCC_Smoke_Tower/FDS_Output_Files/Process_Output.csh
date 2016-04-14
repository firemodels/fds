#!/bin/csh -f
set DDIR=$FDSSMV/Validation/NRCC_Smoke_Tower/Current_Results
set WDIR=$FDSSMV/Validation/NRCC_Smoke_Tower/FDS_Output_Files
cp $DDIR/*devc.csv $WDIR
cp $DDIR/*hrr.csv $WDIR
cp $DDIR/*git.txt $WDIR

