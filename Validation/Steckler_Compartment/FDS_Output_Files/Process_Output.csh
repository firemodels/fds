#!/bin/csh -f
 
set WDIR=$FDSSMV/Validation/Steckler_Compartment/FDS_Output_Files
set DDIR=$FDSSMV/Validation/Steckler_Compartment/Current_Results
cp $DDIR/Steck*line.csv $WDIR
cp $DDIR/*git.txt $WDIR
cp $DDIR/*devc.csv $WDIR

