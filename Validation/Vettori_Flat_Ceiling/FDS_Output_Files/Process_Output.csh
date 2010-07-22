#!/bin/csh -f
setenv SVNROOT ~/FDS-SMV
set DDIR=$SVNROOT/Validation/Vettori_Flat_Ceiling/Current_Results
set WDIR=$SVNROOT/Validation/Vettori_Flat_Ceiling/FDS_Output_Files
cp $DDIR/Vettori*devc.csv $WDIR
cp $DDIR/Vettori*out $WDIR

