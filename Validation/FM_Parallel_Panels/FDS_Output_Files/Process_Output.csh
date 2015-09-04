#!/bin/csh -f
 
set WDIR=$FDSSMV/Validation/FM_Parallel_Panels/FDS_Output_Files
set DDIR=$FDSSMV/Validation/FM_Parallel_Panels/Current_Results
cp $DDIR/FM*line.csv $WDIR
cp $DDIR/*git.txt $WDIR

