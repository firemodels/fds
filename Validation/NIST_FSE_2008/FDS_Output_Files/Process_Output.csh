#!/bin/csh -f
setenv SVNROOT ~/FDS-SMV
set DDIR=$SVNROOT/Validation/NIST_FSE_2008/Current_Results
set WDIR=$SVNROOT/Validation/NIST_FSE_2008/FDS_Output_Files
cp $DDIR/ISO*devc.csv $WDIR
cp $DDIR/ISO*hrr.csv $WDIR

