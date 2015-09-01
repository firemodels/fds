#!/bin/csh -f
 
set WDIR=$FDSSMV/Validation/UMD_Polymers/FDS_Output_Files
set DDIR=$FDSSMV/Validation/UMD_Polymers/Current_Results
cd $DDIR
cp *devc.csv $WDIR
cp *svn.txt $WDIR
cd $WDIR

