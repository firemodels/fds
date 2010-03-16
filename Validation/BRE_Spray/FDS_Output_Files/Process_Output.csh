#!/bin/csh -f
setenv SVNROOT ~/FDS-SMV
set WDIR=$SVNROOT/Validation/BRE_Spray/FDS_Output_Files
set DDIR=$SVNROOT/Validation/BRE_Spray/Current_Results
cd $DDIR
cp BRE_Spray_*_devc.csv $WDIR
cd $WDIR

