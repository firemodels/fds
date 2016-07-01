#!/bin/csh -f
 
set WDIR=$FDSSMV/Validation/UL_NIST_Vents/FDS_Output_Files
set DDIR=$FDSSMV/Validation/UL_NIST_Vents/Current_Results
cp $DDIR/*devc.csv $WDIR
cp $DDIR/*git.txt $WDIR


