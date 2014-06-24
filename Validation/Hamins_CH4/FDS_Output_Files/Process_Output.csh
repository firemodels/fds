#!/bin/csh -f
setenv SVNROOT ~/FDS-SMV
set WDIR=$SVNROOT/Validation/Hamins_CH4/FDS_Output_Files
set DDIR=$SVNROOT/Validation/Hamins_CH4/Current_Results
cp $DDIR/Hamins*devc.csv $WDIR
cp $DDIR/Hamins*svn.txt $WDIR

