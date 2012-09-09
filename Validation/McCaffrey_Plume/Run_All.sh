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

$RUNFDS $INDIR McCaffrey_14_kW 
$RUNFDS $INDIR McCaffrey_22_kW 
$RUNFDS $INDIR McCaffrey_33_kW 
$RUNFDS $INDIR McCaffrey_45_kW 
$RUNFDS $INDIR McCaffrey_57_kW 

$RUNFDS $INDIR McCaffrey_14_kW_coarse
$RUNFDS $INDIR McCaffrey_22_kW_coarse
$RUNFDS $INDIR McCaffrey_33_kW_coarse
$RUNFDS $INDIR McCaffrey_45_kW_coarse
$RUNFDS $INDIR McCaffrey_57_kW_coarse

$RUNFDSMPI 4 $INDIR McCaffrey_14_kW_fine
$RUNFDSMPI 4 $INDIR McCaffrey_22_kW_fine
$RUNFDSMPI 4 $INDIR McCaffrey_33_kW_fine
$RUNFDSMPI 4 $INDIR McCaffrey_45_kW_fine
$RUNFDSMPI 4 $INDIR McCaffrey_57_kW_fine

echo FDS cases submitted
