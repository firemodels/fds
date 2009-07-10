#!/bin/csh -f
setenv SVNROOT ~/FDS-SMV
set DDIR=$SVNROOT/Validation/INERIS_H2/Current_Results
set WDIR=$SVNROOT/Validation/INERIS_H2/FDS_Output_Files
cp $DDIR/INERIS_H2_devc.csv $WDIR

