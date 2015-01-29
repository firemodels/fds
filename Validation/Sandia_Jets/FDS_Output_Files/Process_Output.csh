#!/bin/csh -f
setenv SVNROOT ~/FDS-SMV
set DDIR=$SVNROOT/Validation/Sandia_Jets/Current_Results
set WDIR=$SVNROOT/Validation/Sandia_Jets/FDS_Output_Files
#cp $DDIR/*devc.csv $WDIR
cp $DDIR/*line.csv $WDIR
cp $DDIR/*svn.txt $WDIR

