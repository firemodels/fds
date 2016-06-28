#!/bin/csh -f
 
set WDIR=$FDSSMV/Validation/Hamins_CH4/FDS_Output_Files
set DDIR=$FDSSMV/Validation/Hamins_CH4/Current_Results
cp $DDIR/Hamins*devc.csv $WDIR
cp $DDIR/Hamins*git.txt $WDIR

