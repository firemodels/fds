#!/bin/csh -f
 
set DDIR=$FDSSMV/Validation/CSIRO_Grassland_Fires/Current_Results
set WDIR=$FDSSMV/Validation/CSIRO_Grassland_Fires/FDS_Output_Files
cp $DDIR/*devc.csv $WDIR
cp $DDIR/*git.txt $WDIR

