#!/bin/bash
SVNROOT=~/$1
DIR=$2
cd $SVNROOT/FDS_Compilation/$DIR
rm -f *.o *.mod *.i90
./make_fds.sh
