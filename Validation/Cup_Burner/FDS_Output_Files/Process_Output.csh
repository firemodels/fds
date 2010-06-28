#!/bin/csh -f
setenv SVNROOT ~/FDS-SMV
set PDIR=$SVNROOT/Utilities/Data_Processing
set WDIR=$SVNROOT/Validation/Cup_Burner/FDS_Output_Files
set DDIR=$SVNROOT/Validation/Cup_Burner/Current_Results
cd $DDIR
cp *_devc.csv $WDIR
cd $WDIR
