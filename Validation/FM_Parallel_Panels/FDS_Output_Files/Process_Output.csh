#!/bin/csh -f
setenv SVNROOT ~/FDS-SMV
set WDIR=$SVNROOT/Validation/FM_Parallel_Panels/FDS_Output_Files
set DDIR=$SVNROOT/Validation/FM_Parallel_Panels/Current_Results
cp $DDIR/FM*line.csv $WDIR

