#!/bin/csh -f
setenv SVNROOT ~/FDS-SMV
set DDIR=$SVNROOT/Validation/CSIRO_Grassland_Fires/Current_Results
set WDIR=$SVNROOT/Validation/CSIRO_Grassland_Fires/FDS_Output_Files
cp $DDIR/*devc.csv $WDIR
cp $DDIR/*svn.txt $WDIR

