#!/bin/csh -f
setenv SVNROOT ~/FDS-SMV
set DDIR=$SVNROOT/Validation/Turbulent_Jet/Current_Results
set WDIR=$SVNROOT/Validation/Turbulent_Jet/FDS_Output_Files
cp $DDIR/*line.csv $WDIR
cp $DDIR/*svn.txt $WDIR

