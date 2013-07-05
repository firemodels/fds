#!/bin/csh -f
setenv SVNROOT ~/FDS-SMV
set WDIR=$SVNROOT/Validation/CHRISTIFIRE/FDS_Output_Files
set DDIR=$SVNROOT/Validation/CHRISTIFIRE/Current_Results
cd $DDIR
cp CHRISTIFIRE*devc.csv $WDIR
cp *svn.txt $WDIR
cd $WDIR

