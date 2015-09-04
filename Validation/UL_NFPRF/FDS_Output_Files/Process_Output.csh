#!/bin/csh -f
 
set DDIR=$FDSSMV/Validation/UL_NFPRF/Current_Results
set WDIR=$FDSSMV/Validation/UL_NFPRF/FDS_Output_Files
cp $DDIR/UL*devc.csv $WDIR
cp $DDIR/UL*ctrl.csv $WDIR
cp $DDIR/*git.txt $WDIR

