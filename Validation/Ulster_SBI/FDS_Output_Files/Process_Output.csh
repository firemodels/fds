#!/bin/csh -f
setenv SVNROOT ~/FDS-SMV
set WDIR=$SVNROOT/Validation/Ulster_SBI/FDS_Output_Files
set DDIR=$SVNROOT/Validation/Ulster_SBI/Current_Results
cd $DDIR
cp Ulster*line.csv $WDIR

