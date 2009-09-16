#!/bin/csh -f
setenv SVNROOT ~/FDS-SMV
set DDIR=$SVNROOT/Validation/Decaying_Isotropic_Turbulence/Current_Results
set WDIR=$SVNROOT/Validation/Decaying_Isotropic_Turbulence/FDS_Output_Files
cp $DDIR/*devc.csv $WDIR

