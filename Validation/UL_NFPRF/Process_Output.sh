#!/bin/bash
PDIR=../../../out
CUR=`pwd`
DIR=`basename $CUR`
 
DDIR=Current_Results
WDIR=$PDIR/$DIR
cp $DDIR/UL*devc.csv $WDIR
cp $DDIR/UL*ctrl.csv $WDIR
cp $DDIR/*git.txt $WDIR
