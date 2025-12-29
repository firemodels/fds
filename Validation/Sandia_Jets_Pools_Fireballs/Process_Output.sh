#!/bin/bash
PDIR=../../../out
CUR=`pwd`
DIR=`basename $CUR`
 
DDIR=Current_Results
WDIR=$PDIR/$DIR
cp $DDIR/*devc.csv $WDIR
cp $DDIR/*git.txt $WDIR

