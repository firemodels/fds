#!/bin/csh -f
setenv SVNROOT ~/FDS-SMV
set WDIR=$SVNROOT/Validation/Bryant_Doorway/FDS_Output_Files
set DDIR=$SVNROOT/Validation/Bryant_Doorway/Current_Results
cp $DDIR/Br*line.csv $WDIR


