#!/bin/bash

# This script runs a set of Validation Cases on a Linux machine with a batch queuing system.
# See the file Validation/Common_Run_All.sh for more information.
export SVNROOT=`pwd`/../..
source $SVNROOT/Validation/Common_Run_All.sh

# Run FDS cases
$QFDS $DEBUG $QUEUE -d $INDIR Ranz_Marshall_Table_1_1.fds
$QFDS $DEBUG $QUEUE -d $INDIR Ranz_Marshall_Table_1_2.fds
$QFDS $DEBUG $QUEUE -d $INDIR Ranz_Marshall_Table_1_3.fds
$QFDS $DEBUG $QUEUE -d $INDIR Ranz_Marshall_Table_1_4.fds
$QFDS $DEBUG $QUEUE -d $INDIR Ranz_Marshall_Table_1_5.fds
$QFDS $DEBUG $QUEUE -d $INDIR Ranz_Marshall_Table_1_6.fds
$QFDS $DEBUG $QUEUE -d $INDIR Ranz_Marshall_Table_1_7.fds
$QFDS $DEBUG $QUEUE -d $INDIR Ranz_Marshall_Table_1_8.fds
$QFDS $DEBUG $QUEUE -d $INDIR Ranz_Marshall_Table_1_9.fds
$QFDS $DEBUG $QUEUE -d $INDIR Ranz_Marshall_Table_1_10.fds
$QFDS $DEBUG $QUEUE -d $INDIR Ranz_Marshall_Table_1_11.fds
$QFDS $DEBUG $QUEUE -d $INDIR Ranz_Marshall_Table_1_12.fds
$QFDS $DEBUG $QUEUE -d $INDIR Ranz_Marshall_Table_1_13.fds
$QFDS $DEBUG $QUEUE -d $INDIR Ranz_Marshall_Table_1_14.fds
$QFDS $DEBUG $QUEUE -d $INDIR Ranz_Marshall_Table_1_15.fds
$QFDS $DEBUG $QUEUE -d $INDIR Ranz_Marshall_Table_1_16.fds
$QFDS $DEBUG $QUEUE -d $INDIR Ranz_Marshall_Table_1_17.fds
$QFDS $DEBUG $QUEUE -d $INDIR Ranz_Marshall_Table_1_18.fds
$QFDS $DEBUG $QUEUE -d $INDIR Ranz_Marshall_Table_1_19.fds
$QFDS $DEBUG $QUEUE -d $INDIR Ranz_Marshall_Table_2_1.fds
$QFDS $DEBUG $QUEUE -d $INDIR Ranz_Marshall_Table_2_2.fds
$QFDS $DEBUG $QUEUE -d $INDIR Ranz_Marshall_Table_2_3.fds
$QFDS $DEBUG $QUEUE -d $INDIR Ranz_Marshall_Table_2_4.fds
$QFDS $DEBUG $QUEUE -d $INDIR Ranz_Marshall_Table_2_5.fds
$QFDS $DEBUG $QUEUE -d $INDIR Ranz_Marshall_Table_2_6.fds
$QFDS $DEBUG $QUEUE -d $INDIR Ranz_Marshall_Table_2_7.fds
$QFDS $DEBUG $QUEUE -d $INDIR Ranz_Marshall_Table_2_8.fds
$QFDS $DEBUG $QUEUE -d $INDIR Ranz_Marshall_Table_2_9.fds
$QFDS $DEBUG $QUEUE -d $INDIR Ranz_Marshall_Table_3_1.fds
$QFDS $DEBUG $QUEUE -d $INDIR Ranz_Marshall_Table_3_2.fds
$QFDS $DEBUG $QUEUE -d $INDIR Ranz_Marshall_Table_3_3.fds
$QFDS $DEBUG $QUEUE -d $INDIR Ranz_Marshall_Table_3_4.fds
$QFDS $DEBUG $QUEUE -d $INDIR Ranz_Marshall_Table_3_5.fds
$QFDS $DEBUG $QUEUE -d $INDIR Ranz_Marshall_Table_3_6.fds
$QFDS $DEBUG $QUEUE -d $INDIR Ranz_Marshall_Table_3_7.fds
$QFDS $DEBUG $QUEUE -d $INDIR Ranz_Marshall_Table_3_8.fds
$QFDS $DEBUG $QUEUE -d $INDIR Ranz_Marshall_Table_3_9.fds
$QFDS $DEBUG $QUEUE -d $INDIR Ranz_Marshall_Table_4_1.fds
$QFDS $DEBUG $QUEUE -d $INDIR Ranz_Marshall_Table_4_2.fds
$QFDS $DEBUG $QUEUE -d $INDIR Ranz_Marshall_Table_4_3.fds
$QFDS $DEBUG $QUEUE -d $INDIR Ranz_Marshall_Table_4_4.fds
$QFDS $DEBUG $QUEUE -d $INDIR Ranz_Marshall_Table_4_5.fds
$QFDS $DEBUG $QUEUE -d $INDIR Ranz_Marshall_Table_4_6.fds
$QFDS $DEBUG $QUEUE -d $INDIR Ranz_Marshall_Table_4_7.fds
$QFDS $DEBUG $QUEUE -d $INDIR Ranz_Marshall_Table_4_8.fds
$QFDS $DEBUG $QUEUE -d $INDIR Ranz_Marshall_Table_4_9.fds
$QFDS $DEBUG $QUEUE -d $INDIR Ranz_Marshall_Table_4_10.fds
$QFDS $DEBUG $QUEUE -d $INDIR Ranz_Marshall_Table_4_11.fds
$QFDS $DEBUG $QUEUE -d $INDIR Ranz_Marshall_Table_4_12.fds
$QFDS $DEBUG $QUEUE -d $INDIR Ranz_Marshall_Table_4_13.fds
$QFDS $DEBUG $QUEUE -d $INDIR Ranz_Marshall_Time_Dep.fds

echo FDS cases submitted
