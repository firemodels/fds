#/bin/bash -f

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

$RUNFDSMPI 8 $INDIR UWO_test7_case1_180_32
$RUNFDSMPI 8 $INDIR UWO_test7_case1_180_64
$RUNFDSMPI 8 $INDIR UWO_test7_case1_270_32
$RUNFDSMPI 8 $INDIR UWO_test7_case1_270_64
#$RUNFDSMPI 8 $INDIR UWO_test7_case1_32_spectral
#$RUNFDSMPI 8 $INDIR UWO_test7_case1_64_spectral

echo FDS cases submitted
