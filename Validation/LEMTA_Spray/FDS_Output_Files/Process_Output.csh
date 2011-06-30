#!/bin/csh -f
setenv SVNROOT ~/FDS-SMV
set DDIR=$SVNROOT/Validation/LEMTA_Spray/Current_Results
set WDIR=$SVNROOT/Validation/LEMTA_Spray/FDS_Output_Files
cp $DDIR/LEMT*devc.csv $WDIR
cp $DDIR/*svn.txt $WDIR

