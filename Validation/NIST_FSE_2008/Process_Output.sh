#!/bin/bash
PDIR=../../../out
CUR=`pwd`
DIR=`basename $CUR`
 
DDIR=Current_Results
WDIR=$PDIR/$DIR/FDS_Output_Files
cp $DDIR/ISO*devc.csv $WDIR
cp $DDIR/ISO*hrr.csv $WDIR
cp $DDIR/*git.txt $WDIR

