#!/bin/bash
PDIR=../../../out
CUR=`pwd`
DIR=`basename $CUR`
 
DDIR=Current_Results
WDIR=$PDIR/$DIR
cp $DDIR/Hasenheide*devc.csv $WDIR
cp $DDIR/Hasenheide*git.txt $WDIR
