#!/bin/csh -f
setenv SVNROOT ~/FDS-SMV
set DDIR=$SVNROOT/Validation/SP2009_AST/Current_Results
set WDIR=$SVNROOT/Validation/SP2009_AST/FDS_Output_Files
cp $DDIR/SP2009*devc.csv $WDIR

