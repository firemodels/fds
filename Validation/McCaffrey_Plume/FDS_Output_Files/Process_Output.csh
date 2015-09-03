#!/bin/csh -f
 
set WDIR=$FDSSMV/Validation/McCaffrey_Plume/FDS_Output_Files
set DDIR=$FDSSMV/Validation/McCaffrey_Plume/Current_Results
cd $DDIR
cp McC*line.csv $WDIR
cp *git.txt $WDIR

