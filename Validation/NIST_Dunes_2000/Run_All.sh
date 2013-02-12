#!/bin/bash -f

# This script runs the a set of Validation Cases on a linux machine with
# a batch queuing system

export SVNROOT=`pwd`/../..
export FDSMPI=$SVNROOT/FDS_Compilation/mpi_intel_linux_64/fds_mpi_intel_linux_64
export RUNFDSMPI=$SVNROOT/Utilities/Scripts/runfdsmpi.sh
export BASEDIR=`pwd`
export INDIR=Current_Results
source ~/.bashrc_fds intel64

# uncomment following line to stop all cases
#export STOPFDS=1

/bin/sh -c "cp $BASEDIR/FDS_Input_Files/*.fds $BASEDIR/$INDIR"

$RUNFDSMPI 4 $INDIR NIST_Dunes_2000_SDC02
$RUNFDSMPI 4 $INDIR NIST_Dunes_2000_SDC05
$RUNFDSMPI 4 $INDIR NIST_Dunes_2000_SDC07
$RUNFDSMPI 4 $INDIR NIST_Dunes_2000_SDC10
$RUNFDSMPI 4 $INDIR NIST_Dunes_2000_SDC15
$RUNFDSMPI 4 $INDIR NIST_Dunes_2000_SDC33
$RUNFDSMPI 4 $INDIR NIST_Dunes_2000_SDC35
$RUNFDSMPI 4 $INDIR NIST_Dunes_2000_SDC38
$RUNFDSMPI 4 $INDIR NIST_Dunes_2000_SDC39
