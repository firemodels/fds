#!/bin/csh -f
setenv SVNROOT ~/FDS-SMV
#set DDIR=$SVNROOT/Validation/Utiskul_Compartment/Current_Results
set DDIR=$SVNROOT/Validation/Utiskul_Compartment/Test
set WDIR=$SVNROOT/Validation/Utiskul_Compartment/FDS_Output_Files
cp $DDIR/hu*devc.csv $WDIR
#cp $DDIR/hu*line.csv $WDIR
cp $DDIR/*svn.txt $WDIR

