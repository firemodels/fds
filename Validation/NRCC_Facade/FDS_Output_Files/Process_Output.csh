#!/bin/csh -f
setenv SVNROOT ~/FDS-SMV
set WDIR=$SVNROOT/Validation/NRCC_Facade/FDS_Output_Files
set DDIR=$SVNROOT/Validation/NRCC_Facade/Current_Results
cp $DDIR/NRCC*line.csv $WDIR

