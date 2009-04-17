#!/bin/csh -f
setenv SVNROOT ~/FDS-SMV
set WDIR=$SVNROOT/Validation/VU_Ethanol_Pan_Fire/FDS_Output_Files
set DDIR=$SVNROOT/Validation/VU_Ethanol_Pan_Fire/Current_Results
cd $DDIR
cp VU*hrr.csv $WDIR
cd $WDIR


