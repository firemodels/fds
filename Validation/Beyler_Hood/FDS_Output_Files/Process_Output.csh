#!/bin/csh -f
setenv SVNROOT ~/FDS-SMV
set PDIR=$SVNROOT/Utilities/Data_Processing
set WDIR=$SVNROOT/Validation/Beyler_Hood/FDS_Output_Files
set DDIR=$SVNROOT/Validation/Beyler_Hood/Current_Results
cd $DDIR
$PDIR/Beyler_Hood 
cp Beyler_Hood_FDS.csv $WDIR
cd $WDIR

