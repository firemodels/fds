#!/bin/bash -f

export SVNROOT=`pwd`/../..
export FDS=$SVNROOT/FDS_Compilation/intel_linux_64/fds_intel_linux_64
export RUNFDS=$SVNROOT/Utilities/Scripts/runfds.sh
export FDSMPI=$SVNROOT/FDS_Compilation/mpi_intel_linux_64/fds_mpi_intel_linux_64
export RUNFDSMPI=$SVNROOT/Utilities/Scripts/runfdsmpi.sh
export BASEDIR=`pwd`
export INDIR=Current_Results
source ~/.bashrc_fds intel64

# uncomment following line to stop all cases
# export STOPFDS=1

/bin/sh -c "cp $BASEDIR/FDS_Input_Files/*.fds $BASEDIR/$INDIR"

#$RUNFDS $INDIR Lattimer_20_kW_0_degree_coarse

$RUNFDSMPI 27 $INDIR Lattimer_20_kW_0_degree_coarse
$RUNFDSMPI 27 $INDIR Lattimer_20_kW_0_degree
$RUNFDSMPI 27 $INDIR Lattimer_20_kW_0_degree_fine

echo FDS cases submitted
