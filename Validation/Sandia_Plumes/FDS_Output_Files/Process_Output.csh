#!/bin/csh -f
 
set DDIR=$FDSSMV/Validation/Sandia_Plumes/Current_Results
set WDIR=$FDSSMV/Validation/Sandia_Plumes/FDS_Output_Files
cp $DDIR/Sandia*devc.csv $WDIR
cp $DDIR/Sandia*line.csv $WDIR
cp $DDIR/*git.txt $WDIR

