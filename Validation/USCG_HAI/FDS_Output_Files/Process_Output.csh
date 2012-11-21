#!/bin/csh -f
setenv SVNROOT ~/FDS-SMV
set DDIR=$SVNROOT/Validation/USCG_HAI/Current_Results
set WDIR=$SVNROOT/Validation/USCG_HAI/FDS_Output_Files
cp $DDIR/*devc.csv $WDIR
cp $DDIR/*svn.txt $WDIR

