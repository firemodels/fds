#!/bin/csh -f
setenv SVNROOT ~/FDS-SMV
set WDIR=$SVNROOT/Validation/Restivo_Experiment/FDS_Output_Files
set DDIR=$SVNROOT/Validation/Restivo_Experiment/Current_Results
cd $DDIR
cp Res*line.csv $WDIR

