#!/bin/csh -f
 
set DDIR=$FDSSMV/Validation/Sippola_Aerosol_Deposition/Current_Results
set WDIR=$FDSSMV/Validation/Sippola_Aerosol_Deposition/FDS_Output_Files
cp $DDIR/*devc.csv $WDIR
cp $DDIR/*git.txt $WDIR

