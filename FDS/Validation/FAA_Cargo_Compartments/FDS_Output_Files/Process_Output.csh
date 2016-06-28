#!/bin/csh -f
 
set DDIR=$FDSSMV/Validation/FAA_Cargo_Compartments/Current_Results
set WDIR=$FDSSMV/Validation/FAA_Cargo_Compartments/FDS_Output_Files
cp $DDIR/FAA*devc.csv $WDIR
cp $DDIR/FAA*git.txt $WDIR

