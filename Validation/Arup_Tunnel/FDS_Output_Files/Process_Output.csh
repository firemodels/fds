#!/bin/csh -f
setenv SVNROOT ~/FDS-SMV
set DDIR=$SVNROOT/Validation/Arup_Tunnel/Current_Results
set WDIR=$SVNROOT/Validation/Arup_Tunnel/FDS_Output_Files
cp $DDIR/Arup*devc.csv $WDIR
cp $DDIR/*svn.txt $WDIR

