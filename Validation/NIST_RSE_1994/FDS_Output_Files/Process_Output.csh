#!/bin/csh -f
setenv SVNROOT ~/FDS-SMV
set PDIR=$SVNROOT/Utilities/Data_Processing
set WDIR=$SVNROOT/Validation/NIST_RSE_1994/FDS_Output_Files
set DDIR=$SVNROOT/Validation/NIST_RSE_1994/Current_Results
cd $DDIR
$PDIR/NIST_RSE_1994 
cp NIST_RSE_1994_FDS.csv $WDIR
cd $WDIR

