#!/bin/csh -f
setenv SVNROOT ~/FDS-SMV
set DDIR=$SVNROOT/Validation/LLNL_Enclosure/Current_Results
set WDIR=$SVNROOT/Validation/LLNL_Enclosure/FDS_Output_Files
cp $DDIR/*devc.csv $WDIR

