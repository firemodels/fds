#!/bin/csh -f
 
set WDIR=$FDSSMV/Validation/Waterloo_Methanol/FDS_Output_Files
set DDIR=$FDSSMV/Validation/Waterloo_Methanol/Current_Results
cp $DDIR/*line.csv $WDIR
cp $DDIR/*git.txt $WDIR

