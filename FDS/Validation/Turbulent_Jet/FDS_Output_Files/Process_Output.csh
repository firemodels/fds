#!/bin/csh -f
 
set DDIR=$FDSSMV/Validation/Turbulent_Jet/Current_Results
set WDIR=$FDSSMV/Validation/Turbulent_Jet/FDS_Output_Files
cp $DDIR/*line.csv $WDIR
cp $DDIR/*git.txt $WDIR

