#!/bin/csh -f
 
set DDIR=$FDSSMV/Validation/Wind_Engineering/Current_Results
set WDIR=$FDSSMV/Validation/Wind_Engineering/FDS_Output_Files

cp $DDIR/*line.csv $WDIR
cp $DDIR/*git.txt $WDIR

