#!/bin/bash -f

export SVNROOT=`pwd`/../..
export BASEDIR=`pwd`
export INDIR=Current_Results

source ../../FDS_Compilation/SET_MYFDSENV.sh intel64 run
# to override FDSMPI defined in above script, remove comment
# from  line below and define your own FDSMPI location
#export FDSMPI=override FDSMPI defined in above script

# uncomment following line to stop all cases
#export STOPFDS=1

/bin/sh -c "cp $BASEDIR/FDS_Input_Files/*.fds $BASEDIR/$INDIR"

$RUNFDS $INDIR USCG_HAI_250_kW_Closed_Grinnell
$RUNFDS $INDIR USCG_HAI_500_kW_Closed_Grinnell
$RUNFDS $INDIR USCG_HAI_1000_kW_Closed_Grinnell
$RUNFDS $INDIR USCG_HAI_250_kW_Natural_Grinnell
$RUNFDS $INDIR USCG_HAI_500_kW_Natural_Grinnell
$RUNFDS $INDIR USCG_HAI_1000_kW_Natural_Grinnell
$RUNFDS $INDIR USCG_HAI_250_kW_Forced_Grinnell
$RUNFDS $INDIR USCG_HAI_500_kW_Forced_Grinnell
$RUNFDS $INDIR USCG_HAI_1000_kW_Forced_Grinnell

