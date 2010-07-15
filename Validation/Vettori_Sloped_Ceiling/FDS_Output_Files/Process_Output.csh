#!/bin/csh -f
setenv SVNROOT ~/FDS-SMV
set DDIR=$SVNROOT/Validation/Vettori_Sloped_Ceiling/Current_Results
set WDIR=$SVNROOT/Validation/Vettori_Sloped_Ceiling/FDS_Output_Files
cp $DDIR/Vettori*devc.csv $WDIR

