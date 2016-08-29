#!/bin/bash
 
DDIR=Current_Results
WDIR=FDS_Output_Files
cp $DDIR/UL*devc.csv $WDIR
cp $DDIR/UL*ctrl.csv $WDIR
cp $DDIR/*git.txt $WDIR
