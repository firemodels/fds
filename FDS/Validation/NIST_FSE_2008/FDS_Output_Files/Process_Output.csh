#!/bin/csh -f
 
set DDIR=$FDSSMV/Validation/NIST_FSE_2008/Current_Results
set WDIR=$FDSSMV/Validation/NIST_FSE_2008/FDS_Output_Files
cp $DDIR/ISO*devc.csv $WDIR
cp $DDIR/ISO*hrr.csv $WDIR
cp $DDIR/*git.txt $WDIR

