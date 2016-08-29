#!/bin/bash
 
DDIR=Current_Results
WDIR=FDS_Output_Files
cp $DDIR/*ctrl.csv $WDIR
cp $DDIR/*devc.csv $WDIR
cp $DDIR/*git.txt $WDIR

