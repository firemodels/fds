#!/bin/csh -f
 
set WDIR=$FDSSMV/Validation/NIST_He_2009/FDS_Output_Files
set DDIR=$FDSSMV/Validation/NIST_He_2009/Current_Results
cp $DDIR/*devc.csv $WDIR
cp $DDIR/*git.txt $WDIR

