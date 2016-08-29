#!/bin/bash
 
WDIR=FDS_Output_Files
DDIR=Current_Results
cp $DDIR/*line.csv $WDIR/.
cp $DDIR/*devc.csv $WDIR/.
cp $DDIR/*git.txt $WDIR/.

