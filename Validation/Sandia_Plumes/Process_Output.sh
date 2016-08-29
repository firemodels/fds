#!/bin/bash
 
DDIR=Current_Results
WDIR=FDS_Output_Files
cp $DDIR/Sandia*devc.csv $WDIR
cp $DDIR/Sandia*line.csv $WDIR
cp $DDIR/*git.txt $WDIR

