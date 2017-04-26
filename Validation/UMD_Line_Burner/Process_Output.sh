#!/bin/bash
PDIR=../../../out
CUR=`pwd`
DIR=`basename $CUR`

DDIR=Current_Results
WDIR=$PDIR/$DIR/FDS_Output_Files

cp $DDIR/*git.txt $WDIR

cp $DDIR/methane_dx*line.csv $WDIR

cp $DDIR/methane_XO2*hrr.csv $WDIR
cp $DDIR/propane_XO2*hrr.csv $WDIR

cp $DDIR/methane_XO2*devc.csv $WDIR
cp $DDIR/propane_XO2*devc.csv $WDIR

