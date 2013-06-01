#!/bin/bash

export SVNROOT=`pwd`/../..
export QFDS=/usr/local/bin/qfds.sh
export BASEDIR=`pwd`
export INDIR=Current_Results
# qq="-q fire80s"
qq=

# uncomment following line to stop all cases
# export STOPFDS=1

/bin/sh -c "cp $BASEDIR/FDS_Input_Files/*.fds $BASEDIR/$INDIR"

$QFDS -r $qq -d $INDIR USCG_HAI_250_kW_Closed_Grinnell.fds
$QFDS -r $qq -d $INDIR USCG_HAI_500_kW_Closed_Grinnell.fds
$QFDS -r $qq -d $INDIR USCG_HAI_1000_kW_Closed_Grinnell.fds
$QFDS -r $qq -d $INDIR USCG_HAI_250_kW_Natural_Grinnell.fds
$QFDS -r $qq -d $INDIR USCG_HAI_500_kW_Natural_Grinnell.fds
$QFDS -r $qq -d $INDIR USCG_HAI_1000_kW_Natural_Grinnell.fds
$QFDS -r $qq -d $INDIR USCG_HAI_250_kW_Forced_Grinnell.fds
$QFDS -r $qq -d $INDIR USCG_HAI_500_kW_Forced_Grinnell.fds
$QFDS -r $qq -d $INDIR USCG_HAI_1000_kW_Forced_Grinnell.fds

$QFDS -r $qq -d $INDIR USCG_HAI_250_kW_Closed_Navy.fds
$QFDS -r $qq -d $INDIR USCG_HAI_500_kW_Closed_Navy.fds
$QFDS -r $qq -d $INDIR USCG_HAI_1000_kW_Closed_Navy.fds
$QFDS -r $qq -d $INDIR USCG_HAI_250_kW_Natural_Navy.fds
$QFDS -r $qq -d $INDIR USCG_HAI_500_kW_Natural_Navy.fds
$QFDS -r $qq -d $INDIR USCG_HAI_1000_kW_Natural_Navy.fds
$QFDS -r $qq -d $INDIR USCG_HAI_250_kW_Forced_Navy.fds
$QFDS -r $qq -d $INDIR USCG_HAI_500_kW_Forced_Navy.fds
$QFDS -r $qq -d $INDIR USCG_HAI_1000_kW_Forced_Navy.fds

$QFDS -r $qq -d $INDIR USCG_HAI_250_kW_Closed_Fike.fds
$QFDS -r $qq -d $INDIR USCG_HAI_500_kW_Closed_Fike.fds
$QFDS -r $qq -d $INDIR USCG_HAI_1000_kW_Closed_Fike.fds
$QFDS -r $qq -d $INDIR USCG_HAI_250_kW_Natural_Fike.fds
$QFDS -r $qq -d $INDIR USCG_HAI_500_kW_Natural_Fike.fds
$QFDS -r $qq -d $INDIR USCG_HAI_1000_kW_Natural_Fike.fds
$QFDS -r $qq -d $INDIR USCG_HAI_250_kW_Forced_Fike.fds
$QFDS -r $qq -d $INDIR USCG_HAI_500_kW_Forced_Fike.fds
$QFDS -r $qq -d $INDIR USCG_HAI_1000_kW_Forced_Fike.fds

$QFDS -r $qq -d $INDIR USCG_HAI_250_kW_Closed_Fogtec.fds
$QFDS -r $qq -d $INDIR USCG_HAI_500_kW_Closed_Fogtec.fds
$QFDS -r $qq -d $INDIR USCG_HAI_1000_kW_Closed_Fogtec.fds
$QFDS -r $qq -d $INDIR USCG_HAI_250_kW_Natural_Fogtec.fds
$QFDS -r $qq -d $INDIR USCG_HAI_500_kW_Natural_Fogtec.fds
$QFDS -r $qq -d $INDIR USCG_HAI_1000_kW_Natural_Fogtec.fds
$QFDS -r $qq -d $INDIR USCG_HAI_250_kW_Forced_Fogtec.fds
$QFDS -r $qq -d $INDIR USCG_HAI_500_kW_Forced_Fogtec.fds
$QFDS -r $qq -d $INDIR USCG_HAI_1000_kW_Forced_Fogtec.fds

$QFDS -r $qq -d $INDIR USCG_HAI_WS_29.fds
$QFDS -r $qq -d $INDIR USCG_HAI_WS_45.fds
$QFDS -r $qq -d $INDIR USCG_HAI_WS_50.fds

