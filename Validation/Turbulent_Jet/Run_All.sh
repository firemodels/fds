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

$RUNFDS $INDIR jet_csmag_dx10cm
$RUNFDS $INDIR jet_dsmag_dx10cm
$RUNFDS $INDIR jet_deardorff_dx10cm
$RUNFDS $INDIR jet_vreman_dx10cm

$RUNFDSMPI 49 $INDIR jet_csmag_dx5cm
$RUNFDSMPI 49 $INDIR jet_dsmag_dx5cm
$RUNFDSMPI 49 $INDIR jet_deardorff_dx5cm
$RUNFDSMPI 49 $INDIR jet_vreman_dx5cm

echo FDS cases submitted
