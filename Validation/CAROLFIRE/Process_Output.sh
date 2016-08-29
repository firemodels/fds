#!/bin/bash
PDIR=..
CUR=`pwd`
DIR=`basename $CUR`
 
DDIR=Current_Results
WDIR=$PDIR/$DIR/FDS_Output_Files
cp $DDIR/CAROLFIRE*devc.csv $WDIR
cp $DDIR/*git.txt $WDIR

