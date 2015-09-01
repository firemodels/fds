#!/bin/csh -f
 
set WDIR=$FDSSMV/Validation/CHRISTIFIRE/FDS_Output_Files
set DDIR=$FDSSMV/Validation/CHRISTIFIRE/Current_Results
cd $DDIR
cp CHRISTIFIRE*devc.csv $WDIR
cp *svn.txt $WDIR
cd $WDIR

