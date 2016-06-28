#!/bin/csh -f
 
set WDIR=$FDSSMV/Validation/FAA_Polymers/FDS_Output_Files
set DDIR=$FDSSMV/Validation/FAA_Polymers/Current_Results
cd $DDIR
cp *devc.csv $WDIR
cp *git.txt $WDIR
cd $WDIR

