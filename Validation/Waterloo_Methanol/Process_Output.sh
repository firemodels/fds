#!/bin/bash
PDIR=../../../out
CUR=`pwd`
DIR=`basename $CUR`
 
WDIR=$PDIR/$DIR/FDS_Output_Files
DDIR=Current_Results
cp $DDIR/*line.csv $WDIR
cp $DDIR/Waterloo_Methanol_Predicted*hrr.csv $WDIR
cp $DDIR/*git.txt $WDIR

