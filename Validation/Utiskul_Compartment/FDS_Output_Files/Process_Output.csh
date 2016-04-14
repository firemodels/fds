#!/bin/csh -f
 
#set DDIR=$FDSSMV/Validation/Utiskul_Compartment/Current_Results
set DDIR=$FDSSMV/Validation/Utiskul_Compartment/Test
set WDIR=$FDSSMV/Validation/Utiskul_Compartment/FDS_Output_Files
cp $DDIR/hu*devc.csv $WDIR
#cp $DDIR/hu*line.csv $WDIR
cp $DDIR/*git.txt $WDIR

