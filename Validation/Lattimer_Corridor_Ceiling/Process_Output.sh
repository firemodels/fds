#!/bin/bash
PDIR=../../../out
CUR=`pwd`
DIR=`basename $CUR`

DDIR=Current_Results
WDIR=$PDIR/$DIR

cp $DDIR/*git.txt $WDIR
cp $DDIR/*line.csv $WDIR

