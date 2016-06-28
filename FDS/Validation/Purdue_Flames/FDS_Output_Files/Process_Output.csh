#!/bin/csh -f
 
set DDIR=$FDSSMV/Validation/Purdue_Flames/Current_Results
set WDIR=$FDSSMV/Validation/Purdue_Flames/FDS_Output_Files
cp $DDIR/*line.csv $WDIR
cp $DDIR/*git.txt $WDIR

