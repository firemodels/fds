#!/bin/csh -f
 
set DDIR=$FDSSMV/Validation/Sandia_Jets/Current_Results
set WDIR=$FDSSMV/Validation/Sandia_Jets/FDS_Output_Files
#cp $DDIR/*devc.csv $WDIR
cp $DDIR/*line.csv $WDIR
cp $DDIR/*git.txt $WDIR

