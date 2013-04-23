#!/bin/csh -f
setenv SVNROOT ~/FDS-SMV
set DDIR=$SVNROOT/Validation/VTT_Sprays/Current_Results
set WDIR=$SVNROOT/Validation/VTT_Sprays/FDS_Output_Files
cp $DDIR/*line.csv $WDIR
cp $DDIR/*svn.txt $WDIR

