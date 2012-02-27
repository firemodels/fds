#!/bin/csh -f
setenv SVNROOT ~/FDS-SMV
set WDIR=$SVNROOT/Training/DuPont/FDS_Output_Files
set DDIR=$SVNROOT/Training/DuPont/Current_Results
cp $DDIR/*line.csv $WDIR
cp $DDIR/*devc.csv $WDIR

