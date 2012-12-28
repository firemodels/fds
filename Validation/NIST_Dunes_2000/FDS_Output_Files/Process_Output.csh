#!/bin/csh -f
setenv SVNROOT ~/FDS-SMV
set DDIR=$SVNROOT/Validation/Sippola_Aerosol_Deposition/Current_Results
set WDIR=$SVNROOT/Validation/Sippola_Aerosol_Deposition/FDS_Output_Files
cp $DDIR/*devc.csv $WDIR
cp $DDIR/*svn.txt $WDIR

