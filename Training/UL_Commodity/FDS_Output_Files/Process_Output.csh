#!/bin/csh -f
 
set DDIR=$FDSSMV/Validation/UL_Commodity/Current_Results
set WDIR=$FDSSMV/Validation/UL_Commodity/FDS_Output_Files
cp $DDIR/Group*hrr.csv $WDIR
cp $DDIR/*git.txt $WDIR

