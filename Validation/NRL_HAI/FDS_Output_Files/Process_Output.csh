#!/bin/csh -f
 
set WDIR=$FDSSMV/Validation/NRL_HAI/FDS_Output_Files
set DDIR=$FDSSMV/Validation/NRL_HAI/Current_Results
cd $DDIR
cp NRL*line.csv $WDIR
cp *git.txt $WDIR
cd $WDIR


