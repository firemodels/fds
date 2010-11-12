#!/bin/csh -f
setenv SVNROOT ~/FDS-SMV
set WDIR=$SVNROOT/Validation/McCaffrey_Plume/FDS_Output_Files
set DDIR=$SVNROOT/Validation/McCaffrey_Plume/Current_Results
cd $DDIR
cp McC*line.csv $WDIR

