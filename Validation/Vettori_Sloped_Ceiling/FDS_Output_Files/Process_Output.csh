#!/bin/csh -f
 
set DDIR=$FDSSMV/Validation/Vettori_Sloped_Ceiling/Current_Results
set WDIR=$FDSSMV/Validation/Vettori_Sloped_Ceiling/FDS_Output_Files
cp $DDIR/Vettori*devc.csv $WDIR
cp $DDIR/Vettori*ctrl.csv $WDIR
cp $DDIR/*git.txt $WDIR
