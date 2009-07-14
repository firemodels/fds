#!/bin/csh -f
setenv SVNROOT ~/FDS-SMV
set DDIR=$SVNROOT/Validation/CAROLFIRE/Current_Results
set WDIR=$SVNROOT/Validation/CAROLFIRE/FDS_Output_Files
cp $DDIR/CAROLFIRE*devc.csv $WDIR

