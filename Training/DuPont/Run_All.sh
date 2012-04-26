#!/bin/bash -f

export SVNROOT=`pwd`/../..
export FDS=~/FDS/FDS5/bin/fds5_linux_64
export FDSMPI=$SVNROOT/FDS_Compilation/mpi_intel_linux_64_7031/fds5_mpi_intel_linux_64
export RUNFDS=$SVNROOT/Utilities/Scripts/runfds.sh
export RUNFDSMPI=$SVNROOT/Utilities/Scripts/runfdsmpi.sh
export BASEDIR=`pwd`
export INDIR=Current_Results
source ~/.bashrc_fds intel64

# uncomment following line to stop all cases
# export STOPFDS=1

/bin/sh -c "cp $BASEDIR/FDS_Input_Files/*.fds $BASEDIR/$INDIR"

$RUNFDS $INDIR ASTM_E_648_22
$RUNFDS $INDIR ASTM_E_648_23
$RUNFDS $INDIR ASTM_E_648_24

$RUNFDSMPI 16 $INDIR ASTM_E_648_25

echo FDS cases submitted
