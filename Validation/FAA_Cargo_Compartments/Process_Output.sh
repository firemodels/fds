#!/bin/bash
PDIR=../../../out
CUR=`pwd`
DIR=`basename $CUR`
 
DDIR=Current_Results
WDIR=$PDIR/$DIR
cp $DDIR/FAA*devc.csv $WDIR
cp $DDIR/FAA*git.txt $WDIR

