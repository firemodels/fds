#!/bin/bash

# This script runs a set of Validation Cases on a Linux machine with a batch queuing system.
# See the file Validation/Common_Run_All.sh for more information.
export SVNROOT=`pwd`/../..
source $SVNROOT/Validation/Common_Run_All.sh

$QFDS $DEBUG $QUEUE -d $INDIR USCG_HAI_250_kW_Closed_Grinnell.fds
$QFDS $DEBUG $QUEUE -d $INDIR USCG_HAI_500_kW_Closed_Grinnell.fds
$QFDS $DEBUG $QUEUE -d $INDIR USCG_HAI_1000_kW_Closed_Grinnell.fds
$QFDS $DEBUG $QUEUE -d $INDIR USCG_HAI_250_kW_Natural_Grinnell.fds
$QFDS $DEBUG $QUEUE -d $INDIR USCG_HAI_500_kW_Natural_Grinnell.fds
$QFDS $DEBUG $QUEUE -d $INDIR USCG_HAI_1000_kW_Natural_Grinnell.fds
$QFDS $DEBUG $QUEUE -d $INDIR USCG_HAI_250_kW_Forced_Grinnell.fds
$QFDS $DEBUG $QUEUE -d $INDIR USCG_HAI_500_kW_Forced_Grinnell.fds
$QFDS $DEBUG $QUEUE -d $INDIR USCG_HAI_1000_kW_Forced_Grinnell.fds

$QFDS $DEBUG $QUEUE -d $INDIR USCG_HAI_250_kW_Closed_Navy.fds
$QFDS $DEBUG $QUEUE -d $INDIR USCG_HAI_500_kW_Closed_Navy.fds
$QFDS $DEBUG $QUEUE -d $INDIR USCG_HAI_1000_kW_Closed_Navy.fds
$QFDS $DEBUG $QUEUE -d $INDIR USCG_HAI_250_kW_Natural_Navy.fds
$QFDS $DEBUG $QUEUE -d $INDIR USCG_HAI_500_kW_Natural_Navy.fds
$QFDS $DEBUG $QUEUE -d $INDIR USCG_HAI_1000_kW_Natural_Navy.fds
$QFDS $DEBUG $QUEUE -d $INDIR USCG_HAI_250_kW_Forced_Navy.fds
$QFDS $DEBUG $QUEUE -d $INDIR USCG_HAI_500_kW_Forced_Navy.fds
$QFDS $DEBUG $QUEUE -d $INDIR USCG_HAI_1000_kW_Forced_Navy.fds

$QFDS $DEBUG $QUEUE -d $INDIR USCG_HAI_250_kW_Closed_Fike.fds
$QFDS $DEBUG $QUEUE -d $INDIR USCG_HAI_500_kW_Closed_Fike.fds
$QFDS $DEBUG $QUEUE -d $INDIR USCG_HAI_1000_kW_Closed_Fike.fds
$QFDS $DEBUG $QUEUE -d $INDIR USCG_HAI_250_kW_Natural_Fike.fds
$QFDS $DEBUG $QUEUE -d $INDIR USCG_HAI_500_kW_Natural_Fike.fds
$QFDS $DEBUG $QUEUE -d $INDIR USCG_HAI_1000_kW_Natural_Fike.fds
$QFDS $DEBUG $QUEUE -d $INDIR USCG_HAI_250_kW_Forced_Fike.fds
$QFDS $DEBUG $QUEUE -d $INDIR USCG_HAI_500_kW_Forced_Fike.fds
$QFDS $DEBUG $QUEUE -d $INDIR USCG_HAI_1000_kW_Forced_Fike.fds

$QFDS $DEBUG $QUEUE -d $INDIR USCG_HAI_250_kW_Closed_Fogtec.fds
$QFDS $DEBUG $QUEUE -d $INDIR USCG_HAI_500_kW_Closed_Fogtec.fds
$QFDS $DEBUG $QUEUE -d $INDIR USCG_HAI_1000_kW_Closed_Fogtec.fds
$QFDS $DEBUG $QUEUE -d $INDIR USCG_HAI_250_kW_Natural_Fogtec.fds
$QFDS $DEBUG $QUEUE -d $INDIR USCG_HAI_500_kW_Natural_Fogtec.fds
$QFDS $DEBUG $QUEUE -d $INDIR USCG_HAI_1000_kW_Natural_Fogtec.fds
$QFDS $DEBUG $QUEUE -d $INDIR USCG_HAI_250_kW_Forced_Fogtec.fds
$QFDS $DEBUG $QUEUE -d $INDIR USCG_HAI_500_kW_Forced_Fogtec.fds
$QFDS $DEBUG $QUEUE -d $INDIR USCG_HAI_1000_kW_Forced_Fogtec.fds

