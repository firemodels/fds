#!/bin/bash
PDIR=../../../out
CUR=`pwd`
DIR=`basename $CUR`
 
DDIR=Current_Results
WDIR=$PDIR/$DIR/FDS_Output_Files
cp $DDIR/Vettori*devc.csv $WDIR
cp $DDIR/Vettori*ctrl.csv $WDIR
cp $DDIR/*git.txt $WDIR

