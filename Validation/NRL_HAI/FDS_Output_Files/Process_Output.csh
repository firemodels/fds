#!/bin/csh -f
setenv SVNROOT ~/FDS-SMV
set WDIR=$SVNROOT/Validation/NRL_HAI/FDS_Output_Files
set DDIR=$SVNROOT/Validation/NRL_HAI/Current_Results
cd $DDIR
cp NRL*line.csv $WDIR
cd $WDIR


