#!/bin/csh -f
setenv SVNROOT ~/FDS-SMV
set DDIR=$SVNROOT/Validation/Sandia_Plumes/Current_Results
set WDIR=$SVNROOT/Validation/Sandia_Plumes/FDS_Output_Files
cp $DDIR/Sandia*devc.csv $WDIR
cp $DDIR/Sandia*line.csv $WDIR
cp $DDIR/*svn.txt $WDIR

