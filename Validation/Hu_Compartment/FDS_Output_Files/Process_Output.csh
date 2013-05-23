#!/bin/csh -f
setenv SVNROOT ~/FDS-SMV
set DDIR=$SVNROOT/Validation/Hu_Compartment/Current_Results
set WDIR=$SVNROOT/Validation/Hu_Compartment/FDS_Output_Files
cp $DDIR/hu*devc.csv $WDIR
#cp $DDIR/hu*line.csv $WDIR
cp $DDIR/*svn.txt $WDIR

