#!/bin/csh -f
setenv SVNROOT ~/FDS-SMV
set DDIR=$SVNROOT/Validation/Sandia_He_1m/Current_Results
set WDIR=$SVNROOT/Validation/Sandia_He_1m/FDS_Output_Files
cp $DDIR/Sandia*devc.csv $WDIR

