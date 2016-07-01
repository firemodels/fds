#!/bin/csh -f
 
set PDIR=$FDSSMV/Utilities/Data_Processing
set WDIR=$FDSSMV/Validation/Cup_Burner/FDS_Output_Files
set DDIR=$FDSSMV/Validation/Cup_Burner/Current_Results
cd $DDIR
cp *_devc.csv $WDIR
cp *git.txt $WDIR
cd $WDIR
