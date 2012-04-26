#!/bin/bash

# This script runs the a set of Validation Cases on a linux machine with
# a batch queuing system

export SVNROOT=`pwd`/../..
export FDS=~/FDS/FDS5/bin/fds5_linux_64
export FDSMPI=$SVNROOT/FDS_Compilation/mpi_intel_linux_64_7031/fds5_mpi_intel_linux_64
export RUNFDS=$SVNROOT/Utilities/Scripts/runfds.sh
export RUNFDSMPI=$SVNROOT/Utilities/Scripts/runfdsmpi.sh
export BASEDIR=`pwd`
export INDIR=Current_Results
source ~/.bashrc_fds intel64

# uncomment following line to stop all cases
#export STOPFDS=1

/bin/sh -c "cp $BASEDIR/FDS_Input_Files/*.fds $BASEDIR/$INDIR"

$RUNFDS $INDIR Annulus
$RUNFDS $INDIR Main_Control_Room_Purge
$RUNFDS $INDIR Main_Control_Room_No_Purge
$RUNFDS $INDIR Pump_Room 

$RUNFDSMPI 7 $INDIR Cable_Spreading_Room
$RUNFDSMPI 4 $INDIR Corridor
$RUNFDSMPI 5 $INDIR Switchgear_Room_Cabinet
$RUNFDSMPI 4 $INDIR Switchgear_Room_MCC
$RUNFDSMPI 2 $INDIR Turbine_Building

