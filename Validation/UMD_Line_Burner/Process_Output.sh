#!/bin/bash
PDIR=../../../out
CUR=`pwd`
DIR=`basename $CUR`

DDIR=Current_Results
WDIR=$PDIR/$DIR/FDS_Output_Files

cp $DDIR/*git.txt $WDIR
cp $DDIR/*line.csv $WDIR
cp $DDIR/*hrr.csv $WDIR
cp $DDIR/*devc.csv $WDIR

