#!/bin/csh -f
 
set WDIR=$FDSSMV/Validation/Fleury_Heat_Flux/FDS_Output_Files
set DDIR=$FDSSMV/Validation/Fleury_Heat_Flux/Current_Results
cd $WDIR
cp $DDIR/*line.csv .
cp $DDIR/*git.txt  .


