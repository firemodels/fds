#!/bin/csh -f
 
set DDIR=$FDSSMV/Validation/Arup_Tunnel/Current_Results
set WDIR=$FDSSMV/Validation/Arup_Tunnel/FDS_Output_Files
cp $DDIR/Arup*devc.csv $WDIR
cp $DDIR/*git.txt $WDIR

