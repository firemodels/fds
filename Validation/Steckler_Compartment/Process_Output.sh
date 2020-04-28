#!/bin/bash
PDIR=../../../out
CUR=`pwd`
DIR=`basename $CUR`
 
WDIR=$PDIR/$DIR
DDIR=Current_Results
cp $DDIR/Steck*line.csv $WDIR
cp $DDIR/*git.txt $WDIR
cp $DDIR/*devc.csv $WDIR

