#!/bin/bash
 
DDIR=Current_Results
WDIR=FDS_Output_Files
cp $DDIR/*devc.csv $WDIR
cp $DDIR/*git.txt $WDIR
cp $DDIR/PRS*hrr.csv $WDIR

