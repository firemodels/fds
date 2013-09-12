#!/bin/csh -f
setenv SVNROOT ~/FDS-SMV
set DDIR=$SVNROOT/Validation/NIST_Smoke_Alarms/Current_Results
set WDIR=$SVNROOT/Validation/NIST_Smoke_Alarms/FDS_Output_Files
cp $DDIR/*ctrl.csv $WDIR
cp $DDIR/*devc.csv $WDIR
cp $DDIR/*svn.txt $WDIR

