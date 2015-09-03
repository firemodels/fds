#!/bin/csh -f
 
set DDIR=$FDSSMV/Validation/NIST_NRC/Current_Results
set WDIR=$FDSSMV/Validation/NIST_NRC/FDS_Output_Files
cp $DDIR/*devc.csv $WDIR
cp $DDIR/*git.txt $WDIR

