#!/bin/csh -f
setenv SVNROOT ~/FDS-SMV
set DDIR=$SVNROOT/Validation/SP_AST/Current_Results
set WDIR=$SVNROOT/Validation/SP_AST/FDS_Output_Files
cp $DDIR/SP*devc.csv $WDIR
cp $DDIR/*svn.txt $WDIR

