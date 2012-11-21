#!/bin/csh -f
setenv SVNROOT ~/FDS-SMV
set DDIR=$SVNROOT/Validation/CSTB_Tunnel/Current_Results
set WDIR=$SVNROOT/Validation/CSTB_Tunnel/FDS_Output_Files
cp $DDIR/CSTB*devc.csv $WDIR
cp $DDIR/*svn.txt $WDIR

