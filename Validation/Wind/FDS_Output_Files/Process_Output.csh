#!/bin/csh -f
setenv SVNROOT ~/FDS-SMV
set DDIR=$SVNROOT/Validation/Wind/Current_Results
set WDIR=$SVNROOT/Validation/Wind/FDS_Output_Files
cp $DDIR/*devc.csv $WDIR
cp $DDIR/*line.csv $WDIR
cp $DDIR/*svn.txt $WDIR

