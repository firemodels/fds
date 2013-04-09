#!/bin/bash -f

export SVNROOT=`pwd`/../..
export FDS=$SVNROOT/FDS_Compilation/intel_linux_64/fds_intel_linux_64
export RUNFDS=$SVNROOT/Utilities/Scripts/runfds.sh
export FDSMPI=$SVNROOT/FDS_Compilation/mpi_intel_linux_64/fds_mpi_intel_linux_64
export RUNFDSMPI=$SVNROOT/Utilities/Scripts/runfdsmpi.sh
export BASEDIR=`pwd`
export INDIR=Current_Results

# uncomment following line to stop all cases
# export STOPFDS=1

/bin/sh -c "cp $BASEDIR/FDS_Input_Files/*.fds $BASEDIR/$INDIR"

$RUNFDS $INDIR 7p1_cm_methane_4mm
$RUNFDS $INDIR 7p1_cm_methane_2mm
#$RUNFDSMPI 16 $INDIR 7p1_cm_methane_2mm_16mesh
#$RUNFDSMPI 128 $INDIR 7p1_cm_methane_1mm_128mesh_les
#$RUNFDSMPI 128 $INDIR 7p1_cm_methane_1mm_128mesh_dns

echo FDS cases submitted
