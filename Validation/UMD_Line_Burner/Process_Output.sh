#!/bin/bash
PDIR=../../../out
CUR=`pwd`
DIR=`basename $CUR`
 
DDIR=Current_Results
WDIR=$PDIR/$DIR/FDS_Output_Files
cp $DDIR/*devc.csv $WDIR
cp $DDIR/*line.csv $WDIR
cp $DDIR/*git.txt $WDIR

