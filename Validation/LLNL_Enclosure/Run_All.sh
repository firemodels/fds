#!/bin/bash -f

export SVNROOT=`pwd`/../..
export QFDS=/usr/local/bin/qfds.sh
export BASEDIR=`pwd`
export INDIR=Current_Results
# qq="-q fire80s"
qq=

# uncomment following line to stop all cases
# export STOPFDS=1

/bin/sh -c "cp $BASEDIR/FDS_Input_Files/*.fds $BASEDIR/$INDIR"

$QFDS -r $qq -d $INDIR LLNL_01.fds
$QFDS -r $qq -d $INDIR LLNL_02.fds 
$QFDS -r $qq -d $INDIR LLNL_03.fds
$QFDS -r $qq -d $INDIR LLNL_04.fds 
$QFDS -r $qq -d $INDIR LLNL_05.fds 
$QFDS -r $qq -d $INDIR LLNL_06.fds 
$QFDS -r $qq -d $INDIR LLNL_07.fds 
$QFDS -r $qq -d $INDIR LLNL_08.fds 
$QFDS -r $qq -d $INDIR LLNL_09.fds 
$QFDS -r $qq -d $INDIR LLNL_10.fds 
$QFDS -r $qq -d $INDIR LLNL_11.fds 
$QFDS -r $qq -d $INDIR LLNL_12.fds
$QFDS -r $qq -d $INDIR LLNL_13.fds
$QFDS -r $qq -d $INDIR LLNL_14.fds
$QFDS -r $qq -d $INDIR LLNL_15.fds
$QFDS -r $qq -d $INDIR LLNL_16.fds
$QFDS -r $qq -d $INDIR LLNL_17.fds
$QFDS -r $qq -d $INDIR LLNL_18.fds
$QFDS -r $qq -d $INDIR LLNL_19.fds
$QFDS -r $qq -d $INDIR LLNL_20.fds
$QFDS -r $qq -d $INDIR LLNL_21.fds
$QFDS -r $qq -d $INDIR LLNL_22.fds
$QFDS -r $qq -d $INDIR LLNL_23.fds
$QFDS -r $qq -d $INDIR LLNL_24.fds
$QFDS -r $qq -d $INDIR LLNL_25.fds
$QFDS -r $qq -d $INDIR LLNL_26.fds
$QFDS -r $qq -d $INDIR LLNL_27.fds
$QFDS -r $qq -d $INDIR LLNL_28.fds
$QFDS -r $qq -d $INDIR LLNL_29.fds
$QFDS -r $qq -d $INDIR LLNL_30.fds
$QFDS -r $qq -d $INDIR LLNL_31.fds
$QFDS -r $qq -d $INDIR LLNL_32.fds
$QFDS -r $qq -d $INDIR LLNL_33.fds
$QFDS -r $qq -d $INDIR LLNL_34.fds
$QFDS -r $qq -d $INDIR LLNL_35.fds
$QFDS -r $qq -d $INDIR LLNL_36.fds
$QFDS -r $qq -d $INDIR LLNL_37.fds
$QFDS -r $qq -d $INDIR LLNL_38.fds
$QFDS -r $qq -d $INDIR LLNL_39.fds
$QFDS -r $qq -d $INDIR LLNL_40.fds
$QFDS -r $qq -d $INDIR LLNL_41.fds
$QFDS -r $qq -d $INDIR LLNL_42.fds
$QFDS -r $qq -d $INDIR LLNL_43.fds
$QFDS -r $qq -d $INDIR LLNL_44.fds
$QFDS -r $qq -d $INDIR LLNL_45.fds
$QFDS -r $qq -d $INDIR LLNL_46.fds
$QFDS -r $qq -d $INDIR LLNL_47.fds
$QFDS -r $qq -d $INDIR LLNL_48.fds
$QFDS -r $qq -d $INDIR LLNL_49.fds
$QFDS -r $qq -d $INDIR LLNL_50.fds
$QFDS -r $qq -d $INDIR LLNL_51.fds
$QFDS -r $qq -d $INDIR LLNL_52.fds
$QFDS -r $qq -d $INDIR LLNL_53.fds
$QFDS -r $qq -d $INDIR LLNL_54.fds
$QFDS -r $qq -d $INDIR LLNL_55.fds
$QFDS -r $qq -d $INDIR LLNL_56.fds
$QFDS -r $qq -d $INDIR LLNL_57.fds
$QFDS -r $qq -d $INDIR LLNL_58.fds
$QFDS -r $qq -d $INDIR LLNL_59.fds
$QFDS -r $qq -d $INDIR LLNL_60.fds
$QFDS -r $qq -d $INDIR LLNL_61.fds
$QFDS -r $qq -d $INDIR LLNL_62.fds
$QFDS -r $qq -d $INDIR LLNL_63.fds
$QFDS -r $qq -d $INDIR LLNL_64.fds
 
echo FDS cases submitted
