#!/bin/csh -f
 
set WDIR=$FDSSMV/Validation/Bouchair_Solar_Chimney/FDS_Output_Files
set DDIR=$FDSSMV/Validation/Bouchair_Solar_Chimney/Current_Results
cp $DDIR/*devc.csv $WDIR
cp $DDIR/*git.txt $WDIR


