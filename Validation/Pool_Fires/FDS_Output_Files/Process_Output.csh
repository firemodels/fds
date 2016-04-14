#!/bin/csh -f
 
set DDIR=$FDSSMV/Validation/Pool_Fires/Current_Results
set WDIR=$FDSSMV/Validation/Pool_Fires/FDS_Output_Files
cp $DDIR/*hrr.csv $WDIR
cp $DDIR/*git.txt $WDIR

