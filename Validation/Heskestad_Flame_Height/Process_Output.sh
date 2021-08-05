#!/bin/bash
PDIR=../../../out
CUR=`pwd`
DIR=`basename $CUR`
 
WDIR=$PDIR/$DIR
DDIR=Current_Results
cp $DDIR/Qs*devc.csv $WDIR
cp $DDIR/Qs*line.csv $WDIR
cp $DDIR/Qs*hrr.csv  $WDIR
cp $DDIR/*git.txt    $WDIR
