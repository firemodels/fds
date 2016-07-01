#!/bin/csh -f
 
set WDIR=$FDSSMV/Validation/LLNL_Enclosure/FDS_Output_Files
set DDIR=$FDSSMV/Validation/LLNL_Enclosure/Current_Results
cd $DDIR
cp *devc.csv $WDIR
cp *git.txt $WDIR
cd $WDIR


