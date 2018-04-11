#!/bin/bash

# This script runs a set of Validation Cases on a Linux machine with a batch queuing system.
# See the file Validation/Common_Run_All.sh for more information.
export SVNROOT=`pwd`/../..
source $SVNROOT/Validation/Common_Run_All.sh
nprocs=12

#$QFDS $DEBUG $QUEUE -p 864 -n $nprocs -d $INDIR strong_scaling_test_864.fds
#$QFDS $DEBUG $QUEUE -p 864 -n $nprocs -d $INDIR weak_scaling_test_864.fds

$QFDS $DEBUG $QUEUE -p 432 -n $nprocs -d $INDIR weak_scaling_test_432.fds
$QFDS $DEBUG $QUEUE -p 288 -n $nprocs -d $INDIR weak_scaling_test_288.fds
$QFDS $DEBUG $QUEUE -p 192 -n $nprocs -d $INDIR weak_scaling_test_192.fds
$QFDS $DEBUG $QUEUE -p 128 -n $nprocs -d $INDIR weak_scaling_test_128.fds
$QFDS $DEBUG $QUEUE -p 64  -n $nprocs -d $INDIR weak_scaling_test_064.fds
$QFDS $DEBUG $QUEUE -p 32  -n $nprocs -d $INDIR weak_scaling_test_032.fds
$QFDS $DEBUG $QUEUE -p 16  -n $nprocs -d $INDIR weak_scaling_test_016.fds
$QFDS $DEBUG $QUEUE -p 8   -n $nprocs -d $INDIR weak_scaling_test_008.fds
$QFDS $DEBUG $QUEUE -p 4   -n $nprocs -d $INDIR weak_scaling_test_004.fds
$QFDS $DEBUG $QUEUE -p 2   -n $nprocs -d $INDIR weak_scaling_test_002.fds
$QFDS $DEBUG $QUEUE -p 1   -n $nprocs -d $INDIR weak_scaling_test_001.fds

$QFDS $DEBUG $QUEUE -p 432 -n $nprocs -d $INDIR strong_scaling_test_432.fds
$QFDS $DEBUG $QUEUE -p 288 -n $nprocs -d $INDIR strong_scaling_test_288.fds
$QFDS $DEBUG $QUEUE -p 192 -n $nprocs -d $INDIR strong_scaling_test_192.fds
$QFDS $DEBUG $QUEUE -p 96  -n $nprocs -d $INDIR strong_scaling_test_096.fds
$QFDS $DEBUG $QUEUE -p 64  -n $nprocs -d $INDIR strong_scaling_test_064.fds
$QFDS $DEBUG $QUEUE -p 32  -n $nprocs -d $INDIR strong_scaling_test_032.fds
$QFDS $DEBUG $QUEUE -p 8   -n $nprocs -d $INDIR strong_scaling_test_008.fds
$QFDS $DEBUG $QUEUE -p 1   -n $nprocs -d $INDIR strong_scaling_test_001.fds

# run GLMAT pressure solver cases

$QFDS $DEBUG $QUEUE -p 128 -n $nprocs -d $INDIR weak_scaling_glmat_128.fds
$QFDS $DEBUG $QUEUE -p 64  -n $nprocs -d $INDIR weak_scaling_glmat_064.fds
$QFDS $DEBUG $QUEUE -p 32  -n $nprocs -d $INDIR weak_scaling_glmat_032.fds
$QFDS $DEBUG $QUEUE -p 16  -n $nprocs -d $INDIR weak_scaling_glmat_016.fds
$QFDS $DEBUG $QUEUE -p 8   -n $nprocs -d $INDIR weak_scaling_glmat_008.fds
$QFDS $DEBUG $QUEUE -p 4   -n $nprocs -d $INDIR weak_scaling_glmat_004.fds
$QFDS $DEBUG $QUEUE -p 2   -n $nprocs -d $INDIR weak_scaling_glmat_002.fds
$QFDS $DEBUG $QUEUE -p 1   -n $nprocs -d $INDIR weak_scaling_glmat_001.fds

$QFDS $DEBUG $QUEUE -p 192 -n $nprocs -d $INDIR strong_scaling_glmat_192.fds
$QFDS $DEBUG $QUEUE -p 96  -n $nprocs -d $INDIR strong_scaling_glmat_096.fds
$QFDS $DEBUG $QUEUE -p 64  -n $nprocs -d $INDIR strong_scaling_glmat_064.fds
$QFDS $DEBUG $QUEUE -p 32  -n $nprocs -d $INDIR strong_scaling_glmat_032.fds
$QFDS $DEBUG $QUEUE -p 8   -n $nprocs -d $INDIR strong_scaling_glmat_008.fds
$QFDS $DEBUG $QUEUE -p 1   -n $nprocs -d $INDIR strong_scaling_glmat_001.fds

# run ScaRC pressure solver cases

$QFDS $DEBUG $QUEUE -p 432 -n $nprocs -d $INDIR weak_scaling_scarc_432.fds
$QFDS $DEBUG $QUEUE -p 288 -n $nprocs -d $INDIR weak_scaling_scarc_288.fds
$QFDS $DEBUG $QUEUE -p 192 -n $nprocs -d $INDIR weak_scaling_scarc_192.fds
$QFDS $DEBUG $QUEUE -p 128 -n $nprocs -d $INDIR weak_scaling_scarc_128.fds
$QFDS $DEBUG $QUEUE -p 64  -n $nprocs -d $INDIR weak_scaling_scarc_064.fds
$QFDS $DEBUG $QUEUE -p 32  -n $nprocs -d $INDIR weak_scaling_scarc_032.fds
$QFDS $DEBUG $QUEUE -p 16  -n $nprocs -d $INDIR weak_scaling_scarc_016.fds
$QFDS $DEBUG $QUEUE -p 8   -n $nprocs -d $INDIR weak_scaling_scarc_008.fds
$QFDS $DEBUG $QUEUE -p 4   -n $nprocs -d $INDIR weak_scaling_scarc_004.fds
$QFDS $DEBUG $QUEUE -p 2   -n $nprocs -d $INDIR weak_scaling_scarc_002.fds
$QFDS $DEBUG $QUEUE -p 1   -n $nprocs -d $INDIR weak_scaling_scarc_001.fds

$QFDS $DEBUG $QUEUE -p 432 -n $nprocs -d $INDIR strong_scaling_scarc_432.fds
$QFDS $DEBUG $QUEUE -p 288 -n $nprocs -d $INDIR strong_scaling_scarc_288.fds
$QFDS $DEBUG $QUEUE -p 192 -n $nprocs -d $INDIR strong_scaling_scarc_192.fds
$QFDS $DEBUG $QUEUE -p 96  -n $nprocs -d $INDIR strong_scaling_scarc_096.fds
$QFDS $DEBUG $QUEUE -p 64  -n $nprocs -d $INDIR strong_scaling_scarc_064.fds
$QFDS $DEBUG $QUEUE -p 32  -n $nprocs -d $INDIR strong_scaling_scarc_032.fds
$QFDS $DEBUG $QUEUE -p 8   -n $nprocs -d $INDIR strong_scaling_scarc_008.fds
$QFDS $DEBUG $QUEUE -p 1   -n $nprocs -d $INDIR strong_scaling_scarc_001.fds
