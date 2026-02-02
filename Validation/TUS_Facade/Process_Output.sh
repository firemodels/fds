#!/bin/bash
PDIR=../../../out
CUR=`pwd`
DIR=`basename $CUR`
 
WDIR=$PDIR/$DIR
DDIR=Current_Results
cp $DDIR/*line.csv $WDIR
cp $DDIR/*devc.csv $WDIR
cp $DDIR/*hrr.csv $WDIR
cp $DDIR/*git.txt $WDIR
