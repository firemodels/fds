#!/bin/csh -f
setenv SVNROOT ~/FDS-SMV
set DDIR=$SVNROOT/Validation/FAA_Cargo_Compartments/Current_Results
set WDIR=$SVNROOT/Validation/FAA_Cargo_Compartments/FDS_Output_Files
cp $DDIR/FAA*devc.csv $WDIR
cp $DDIR/FAA*svn.txt $WDIR

