#!/bin/bash

# This script runs a set of Validation Cases on a Linux machine with a batch queuing system.
# See the file Validation/Common_Run_All.sh for more information.
export SVNROOT=`pwd`/../..
source $SVNROOT/Validation/Common_Run_All.sh

$QFDS $DEBUG $QUEUE -t -o 1 -d $INDIR openmp_test64a.fds
$QFDS $DEBUG $QUEUE -t -o 2 -d $INDIR openmp_test64b.fds
$QFDS $DEBUG $QUEUE -t -o 3 -d $INDIR openmp_test64c.fds
$QFDS $DEBUG $QUEUE -t -o 4 -d $INDIR openmp_test64d.fds
$QFDS $DEBUG $QUEUE -t -o 5 -d $INDIR openmp_test64e.fds
$QFDS $DEBUG $QUEUE -t -o 6 -d $INDIR openmp_test64f.fds
$QFDS $DEBUG $QUEUE -t -o 7 -d $INDIR openmp_test64g.fds
$QFDS $DEBUG $QUEUE -t -o 8 -d $INDIR openmp_test64h.fds
$QFDS $DEBUG $QUEUE -t -o 1 -d $INDIR openmp_test128a.fds
$QFDS $DEBUG $QUEUE -t -o 2 -d $INDIR openmp_test128b.fds
$QFDS $DEBUG $QUEUE -t -o 3 -d $INDIR openmp_test128c.fds
$QFDS $DEBUG $QUEUE -t -o 4 -d $INDIR openmp_test128d.fds
$QFDS $DEBUG $QUEUE -t -o 5 -d $INDIR openmp_test128e.fds
$QFDS $DEBUG $QUEUE -t -o 6 -d $INDIR openmp_test128f.fds
$QFDS $DEBUG $QUEUE -t -o 7 -d $INDIR openmp_test128g.fds
$QFDS $DEBUG $QUEUE -t -o 8 -d $INDIR openmp_test128h.fds
