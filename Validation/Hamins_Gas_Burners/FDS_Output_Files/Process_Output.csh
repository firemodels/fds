#!/bin/csh -f
 
set WDIR=$FDSSMV/Validation/Hamins_Gas_Burners/FDS_Output_Files
set DDIR=$FDSSMV/Validation/Hamins_Gas_Burners/Current_Results
cp $DDIR/*line.csv $WDIR
cp $DDIR/*git.txt $WDIR

