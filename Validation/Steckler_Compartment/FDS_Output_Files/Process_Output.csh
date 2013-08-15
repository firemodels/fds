#!/bin/csh -f
setenv SVNROOT ~/FDS-SMV
set WDIR=$SVNROOT/Validation/Steckler_Compartment/FDS_Output_Files
set DDIR=$SVNROOT/Validation/Steckler_Compartment/Current_Results
cp $DDIR/Steck*line.csv $WDIR
cp $DDIR/*svn.txt $WDIR
cp $DDIR/*devc.csv $WDIR

