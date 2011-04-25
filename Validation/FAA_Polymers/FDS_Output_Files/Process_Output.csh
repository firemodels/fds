#!/bin/csh -f
setenv SVNROOT ~/FDS-SMV
set WDIR=$SVNROOT/Validation/FAA_Polymers/FDS_Output_Files
set DDIR=$SVNROOT/Validation/FAA_Polymers/Current_Results
cd $DDIR
cp FAA*devc.csv $WDIR
cd $WDIR

