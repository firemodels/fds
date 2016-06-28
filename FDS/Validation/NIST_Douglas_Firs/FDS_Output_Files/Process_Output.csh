#!/bin/csh -f
 
set DDIR=$FDSSMV/Validation/NIST_Douglas_Firs/Current_Results
set WDIR=$FDSSMV/Validation/NIST_Douglas_Firs/FDS_Output_Files
cp $DDIR/*hrr.csv $WDIR
cp $DDIR/*git.txt $WDIR

