#!/bin/csh -f
 
set WDIR=$FDSSMV/Validation/MPI_Scaling_Tests/FDS_Output_Files
set DDIR=$FDSSMV/Validation/MPI_Scaling_Tests/Current_Results
cp $DDIR/*cpu.csv $WDIR
cp $DDIR/*git.txt $WDIR


