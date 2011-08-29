#!/bin/csh -f
setenv SVNROOT ~/FDS-SMV
set WDIR=$SVNROOT/Validation/Fleury_Heat_Flux/FDS_Output_Files
set DDIR=$SVNROOT/Validation/Fleury_Heat_Flux/Current_Results
cd $WDIR
cp $DDIR/*line.csv .
cp $DDIR/*svn.txt  .


