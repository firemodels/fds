#!/bin/csh -f
 
set DDIR=$FDSSMV/Validation/VTT_Sprays/Current_Results
set WDIR=$FDSSMV/Validation/VTT_Sprays/FDS_Output_Files
cp $DDIR/*line.csv $WDIR
cp $DDIR/*git.txt $WDIR

