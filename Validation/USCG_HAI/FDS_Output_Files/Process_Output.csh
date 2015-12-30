#!/bin/csh -f
 
set DDIR=$FDSSMV/Validation/USCG_HAI/Current_Results
set WDIR=$FDSSMV/Validation/USCG_HAI/FDS_Output_Files
cp $DDIR/*devc.csv $WDIR
cp $DDIR/*git.txt $WDIR

