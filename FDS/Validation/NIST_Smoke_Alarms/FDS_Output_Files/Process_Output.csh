#!/bin/csh -f
 
set DDIR=$FDSSMV/Validation/NIST_Smoke_Alarms/Current_Results
set WDIR=$FDSSMV/Validation/NIST_Smoke_Alarms/FDS_Output_Files
cp $DDIR/*ctrl.csv $WDIR
cp $DDIR/*devc.csv $WDIR
cp $DDIR/*git.txt $WDIR

