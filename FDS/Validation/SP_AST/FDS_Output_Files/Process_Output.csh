#!/bin/csh -f
 
set DDIR=$FDSSMV/Validation/SP_AST/Current_Results
set WDIR=$FDSSMV/Validation/SP_AST/FDS_Output_Files
cp $DDIR/*devc.csv $WDIR
cp $DDIR/*git.txt $WDIR

