#!/bin/bash
PDIR=../../../out
CUR=`pwd`
DIR=`basename $CUR`
 
WDIR=$PDIR/$DIR
DDIR=Current_Results
cp $DDIR/BRE_Spray_*_devc.csv $WDIR
cp $DDIR/*git.txt $WDIR
