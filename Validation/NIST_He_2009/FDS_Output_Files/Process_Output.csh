#!/bin/csh -f
setenv SVNROOT ~/FDS-SMV
set DDIR=$SVNROOT/Validation/NIST_He_2009/Current_Results
set WDIR=$SVNROOT/Validation/NIST_He_2009/FDS_Output_Files
cp $DDIR/NIST*devc.csv $WDIR

