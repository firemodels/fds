#!/bin/csh -f
setenv SVNROOT ~/FDS-SMV
set DDIR=$SVNROOT/Validation/NIST_Dunes_2000/Current_Results
set WDIR=$SVNROOT/Validation/NIST_Dunes_2000/FDS_Output_Files
cp $DDIR/*ctrl.csv $WDIR
cp $DDIR/*devc.csv $WDIR
cp $DDIR/*svn.txt $WDIR

