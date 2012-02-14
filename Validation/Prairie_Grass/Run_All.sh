#!/bin/bash -f

export SVNROOT=`pwd`/../..
export FDSMPI=$SVNROOT/FDS_Compilation/mpi_intel_linux_64/fds_mpi_intel_linux_64
export RUNFDSMPI=$SVNROOT/Utilities/Scripts/runfdsmpi.sh
export BASEDIR=`pwd`
export INDIR=Current_Results
source ~/.bashrc_fds intel64

# uncomment following line to stop all cases
#export STOPFDS=1

/bin/sh -c "cp $BASEDIR/FDS_Input_Files/*.fds $BASEDIR/$INDIR"

$RUNFDSMPI 8 $INDIR Prairie_Grass_Case_21
$RUNFDSMPI 8 $INDIR Prairie_Grass_Case_22
$RUNFDSMPI 8 $INDIR Prairie_Grass_Case_23
$RUNFDSMPI 8 $INDIR Prairie_Grass_Case_24
$RUNFDSMPI 8 $INDIR Prairie_Grass_Case_42
$RUNFDSMPI 8 $INDIR Prairie_Grass_Case_46
$RUNFDSMPI 8 $INDIR Prairie_Grass_Case_55
$RUNFDSMPI 8 $INDIR Prairie_Grass_Case_57 

echo FDS cases submitted
