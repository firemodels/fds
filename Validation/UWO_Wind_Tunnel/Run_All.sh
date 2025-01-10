#/bin/bash

# This script runs a set of Validation Cases on a Linux machine with a batch queuing system.
# See the file Validation/Common_Run_All.sh for more information.
export SVNROOT=`pwd`/../..
source $SVNROOT/Validation/Common_Run_All.sh

$QFDS $DEBUG $QUEUE -p 8     -d $INDIR UWO_Test7_Case1_180_5p00mm.fds
$QFDS $DEBUG $QUEUE -p 64    -d $INDIR UWO_Test7_Case1_180_2p50mm.fds
$QFDS $DEBUG $QUEUE -p 120   -d $INDIR UWO_Test7_Case1_180_1p25mm.fds
$QFDS $DEBUG $QUEUE -p 8     -d $INDIR UWO_Test7_Case1_270_5p00mm.fds
$QFDS $DEBUG $QUEUE -p 64    -d $INDIR UWO_Test7_Case1_270_2p50mm.fds
$QFDS $DEBUG $QUEUE -p 120   -d $INDIR UWO_Test7_Case1_270_1p25mm.fds
$QFDS $DEBUG $QUEUE -p 8     -d $INDIR UWO_SS21_Test6_40ft_0_20mm.fds
$QFDS $DEBUG $QUEUE -p 64    -d $INDIR UWO_SS21_Test6_40ft_0_10mm.fds
$QFDS $DEBUG $QUEUE -p 128   -d $INDIR UWO_SS21_Test6_40ft_0_5mm.fds
$QFDS $DEBUG $QUEUE -p 8     -d $INDIR UWO_SS21_Test6_40ft_45_20mm.fds
$QFDS $DEBUG $QUEUE -p 64    -d $INDIR UWO_SS21_Test6_40ft_45_10mm.fds
$QFDS $DEBUG $QUEUE -p 128   -d $INDIR UWO_SS21_Test6_40ft_45_5mm.fds

echo FDS cases submitted
