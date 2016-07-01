#!/bin/csh -f
 
set PDIR=$FDSSMV/Utilities/Data_Processing
set WDIR=$FDSSMV/Validation/Harrison_Spill_Plumes/FDS_Output_Files
set DDIR=$FDSSMV/Validation/Harrison_Spill_Plumes/Current_Results
cd $DDIR
cp *devc.csv  $WDIR
cp *git.txt   $WDIR
cd $WDIR


