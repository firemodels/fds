#!/bin/bash

# This script runs a set of Validation Cases on a Linux machine with a batch queuing system.
# See the file Validation/Common_Run_All.sh for more information.
export SVNROOT=`pwd`/../..
source $SVNROOT/Validation/Common_Run_All.sh

$QFDS $DEBUG -p 36  $QUEUE -d $INDIR NIST_Backdraft_Methane_25kW_360s_low_coarse.fds
$QFDS $DEBUG -p 36  $QUEUE -d $INDIR NIST_Backdraft_Methane_25kW_390s_low_coarse.fds
$QFDS $DEBUG -p 36  $QUEUE -d $INDIR NIST_Backdraft_Methane_25kW_450s_low_coarse.fds
$QFDS $DEBUG -p 36  $QUEUE -d $INDIR NIST_Backdraft_Methane_37p5kW_225s_low_coarse.fds
$QFDS $DEBUG -p 36  $QUEUE -d $INDIR NIST_Backdraft_Methane_37p5kW_285s_low_coarse.fds
$QFDS $DEBUG -p 36  $QUEUE -d $INDIR NIST_Backdraft_Propane_25kW_210s_low_coarse.fds
$QFDS $DEBUG -p 36  $QUEUE -d $INDIR NIST_Backdraft_Propane_25kW_285s_low_coarse.fds
$QFDS $DEBUG -p 36  $QUEUE -d $INDIR NIST_Backdraft_Propane_16p7kW_255s_low_coarse.fds
$QFDS $DEBUG -p 36  $QUEUE -d $INDIR NIST_Backdraft_Propane_16p7kW_270s_low_coarse.fds
$QFDS $DEBUG -p 36  $QUEUE -d $INDIR NIST_Backdraft_Propane_16p7kW_315s_low_coarse.fds
$QFDS $DEBUG -p 36  $QUEUE -d $INDIR NIST_Backdraft_Methane_25kW_360s_mid_coarse.fds
$QFDS $DEBUG -p 36  $QUEUE -d $INDIR NIST_Backdraft_Methane_25kW_390s_mid_coarse.fds
$QFDS $DEBUG -p 36  $QUEUE -d $INDIR NIST_Backdraft_Methane_25kW_450s_mid_coarse.fds
$QFDS $DEBUG -p 36  $QUEUE -d $INDIR NIST_Backdraft_Methane_37p5kW_225s_mid_coarse.fds
$QFDS $DEBUG -p 36  $QUEUE -d $INDIR NIST_Backdraft_Methane_37p5kW_285s_mid_coarse.fds
$QFDS $DEBUG -p 36  $QUEUE -d $INDIR NIST_Backdraft_Propane_25kW_210s_mid_coarse.fds
$QFDS $DEBUG -p 36  $QUEUE -d $INDIR NIST_Backdraft_Propane_25kW_285s_mid_coarse.fds
$QFDS $DEBUG -p 36  $QUEUE -d $INDIR NIST_Backdraft_Propane_16p7kW_255s_mid_coarse.fds
$QFDS $DEBUG -p 36  $QUEUE -d $INDIR NIST_Backdraft_Propane_16p7kW_270s_mid_coarse.fds
$QFDS $DEBUG -p 36  $QUEUE -d $INDIR NIST_Backdraft_Propane_16p7kW_315s_mid_coarse.fds
$QFDS $DEBUG -p 36  $QUEUE -d $INDIR NIST_Backdraft_Methane_25kW_360s_low_medium.fds
$QFDS $DEBUG -p 36  $QUEUE -d $INDIR NIST_Backdraft_Methane_25kW_390s_low_medium.fds
$QFDS $DEBUG -p 36  $QUEUE -d $INDIR NIST_Backdraft_Methane_25kW_450s_low_medium.fds
$QFDS $DEBUG -p 36  $QUEUE -d $INDIR NIST_Backdraft_Methane_37p5kW_225s_low_medium.fds
$QFDS $DEBUG -p 36  $QUEUE -d $INDIR NIST_Backdraft_Methane_37p5kW_285s_low_medium.fds
$QFDS $DEBUG -p 36  $QUEUE -d $INDIR NIST_Backdraft_Propane_25kW_210s_low_medium.fds
$QFDS $DEBUG -p 36  $QUEUE -d $INDIR NIST_Backdraft_Propane_25kW_285s_low_medium.fds
$QFDS $DEBUG -p 36  $QUEUE -d $INDIR NIST_Backdraft_Propane_16p7kW_255s_low_medium.fds
$QFDS $DEBUG -p 36  $QUEUE -d $INDIR NIST_Backdraft_Propane_16p7kW_270s_low_medium.fds
$QFDS $DEBUG -p 36  $QUEUE -d $INDIR NIST_Backdraft_Propane_16p7kW_315s_low_medium.fds
$QFDS $DEBUG -p 36  $QUEUE -d $INDIR NIST_Backdraft_Methane_25kW_360s_mid_medium.fds
$QFDS $DEBUG -p 36  $QUEUE -d $INDIR NIST_Backdraft_Methane_25kW_390s_mid_medium.fds
$QFDS $DEBUG -p 36  $QUEUE -d $INDIR NIST_Backdraft_Methane_25kW_450s_mid_medium.fds
$QFDS $DEBUG -p 36  $QUEUE -d $INDIR NIST_Backdraft_Methane_37p5kW_225s_mid_medium.fds
$QFDS $DEBUG -p 36  $QUEUE -d $INDIR NIST_Backdraft_Methane_37p5kW_285s_mid_medium.fds
$QFDS $DEBUG -p 36  $QUEUE -d $INDIR NIST_Backdraft_Propane_25kW_210s_mid_medium.fds
$QFDS $DEBUG -p 36  $QUEUE -d $INDIR NIST_Backdraft_Propane_25kW_285s_mid_medium.fds
$QFDS $DEBUG -p 36  $QUEUE -d $INDIR NIST_Backdraft_Propane_16p7kW_255s_mid_medium.fds
$QFDS $DEBUG -p 36  $QUEUE -d $INDIR NIST_Backdraft_Propane_16p7kW_270s_mid_medium.fds
$QFDS $DEBUG -p 36  $QUEUE -d $INDIR NIST_Backdraft_Propane_16p7kW_315s_mid_medium.fds
$QFDS $DEBUG -p 36  $QUEUE -d $INDIR NIST_Backdraft_Methane_25kW_360s_low_fine.fds
$QFDS $DEBUG -p 36  $QUEUE -d $INDIR NIST_Backdraft_Methane_25kW_390s_low_fine.fds
$QFDS $DEBUG -p 36  $QUEUE -d $INDIR NIST_Backdraft_Methane_25kW_450s_low_fine.fds
$QFDS $DEBUG -p 36  $QUEUE -d $INDIR NIST_Backdraft_Methane_37p5kW_225s_low_fine.fds
$QFDS $DEBUG -p 36  $QUEUE -d $INDIR NIST_Backdraft_Methane_37p5kW_285s_low_fine.fds
$QFDS $DEBUG -p 36  $QUEUE -d $INDIR NIST_Backdraft_Propane_25kW_210s_low_fine.fds
$QFDS $DEBUG -p 36  $QUEUE -d $INDIR NIST_Backdraft_Propane_25kW_285s_low_fine.fds
$QFDS $DEBUG -p 36  $QUEUE -d $INDIR NIST_Backdraft_Propane_16p7kW_255s_low_fine.fds
$QFDS $DEBUG -p 36  $QUEUE -d $INDIR NIST_Backdraft_Propane_16p7kW_270s_low_fine.fds
$QFDS $DEBUG -p 36  $QUEUE -d $INDIR NIST_Backdraft_Propane_16p7kW_315s_low_fine.fds
$QFDS $DEBUG -p 36  $QUEUE -d $INDIR NIST_Backdraft_Methane_25kW_360s_mid_fine.fds
$QFDS $DEBUG -p 36  $QUEUE -d $INDIR NIST_Backdraft_Methane_25kW_390s_mid_fine.fds
$QFDS $DEBUG -p 36  $QUEUE -d $INDIR NIST_Backdraft_Methane_25kW_450s_mid_fine.fds
$QFDS $DEBUG -p 36  $QUEUE -d $INDIR NIST_Backdraft_Methane_37p5kW_225s_mid_fine.fds
$QFDS $DEBUG -p 36  $QUEUE -d $INDIR NIST_Backdraft_Methane_37p5kW_285s_mid_fine.fds
$QFDS $DEBUG -p 36  $QUEUE -d $INDIR NIST_Backdraft_Propane_25kW_210s_mid_fine.fds
$QFDS $DEBUG -p 36  $QUEUE -d $INDIR NIST_Backdraft_Propane_25kW_285s_mid_fine.fds
$QFDS $DEBUG -p 36  $QUEUE -d $INDIR NIST_Backdraft_Propane_16p7kW_255s_mid_fine.fds
$QFDS $DEBUG -p 36  $QUEUE -d $INDIR NIST_Backdraft_Propane_16p7kW_270s_mid_fine.fds
$QFDS $DEBUG -p 36  $QUEUE -d $INDIR NIST_Backdraft_Propane_16p7kW_315s_mid_fine.fds

echo FDS cases submitted
