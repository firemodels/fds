#!/bin/bash
 
DDIR=Current_Results
WDIR=FDS_Output_Files
cp $DDIR/Vettori*devc.csv $WDIR
cp $DDIR/Vettori*ctrl.csv $WDIR
cp $DDIR/*git.txt $WDIR
