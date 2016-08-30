#!/bin/bash
PDIR=..
CUR=`pwd`
DIR=`basename $CUR`
 
DDIR=Current_Results
WDIR=$PDIR/$DIR/FDS_Output_Files
cp $DDIR/FAA*devc.csv $WDIR
cp $DDIR/FAA*git.txt $WDIR

