#!/bin/csh -f
 
set WDIR=$FDSSMV/Validation/BRE_Spray/FDS_Output_Files
set DDIR=$FDSSMV/Validation/BRE_Spray/Current_Results
cd $DDIR
cp BRE_Spray_*_devc.csv $WDIR
cp *git.txt $WDIR
cd $WDIR

