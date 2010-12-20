#!/bin/csh -f
setenv SVNROOT ~/FDS-SMV
set DDIR=$SVNROOT/Validation/LLNL_Enclosures/Current_Results
set WDIR=$SVNROOT/Validation/LLNL_Enclosures/FDS_Output_Files
cp $DDIR/*devc.csv $WDIR

