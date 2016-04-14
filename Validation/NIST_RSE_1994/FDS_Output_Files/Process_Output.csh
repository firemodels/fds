#!/bin/csh -f
 
set WDIR=$FDSSMV/Validation/NIST_RSE_1994/FDS_Output_Files
set DDIR=$FDSSMV/Validation/NIST_RSE_1994/Current_Results
cp $DDIR/*devc.csv $WDIR
cp $DDIR/*git.txt $WDIR

