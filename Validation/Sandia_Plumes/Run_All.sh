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

$RUNFDSMPI 16 $INDIR Sandia_He_1m_dx6cm
$RUNFDSMPI 16 $INDIR Sandia_He_1m_dx3cm
$RUNFDSMPI 16 $INDIR Sandia_He_1m_dx1p5cm
$RUNFDSMPI 16 $INDIR Sandia_CH4_1m_Test14_dx6cm
$RUNFDSMPI 16 $INDIR Sandia_CH4_1m_Test14_dx3cm
$RUNFDSMPI 16 $INDIR Sandia_CH4_1m_Test14_dx1p5cm
$RUNFDSMPI 16 $INDIR Sandia_CH4_1m_Test17_dx6cm
$RUNFDSMPI 16 $INDIR Sandia_CH4_1m_Test17_dx3cm
$RUNFDSMPI 16 $INDIR Sandia_CH4_1m_Test17_dx1p5cm
$RUNFDSMPI 16 $INDIR Sandia_CH4_1m_Test24_dx6cm
$RUNFDSMPI 16 $INDIR Sandia_CH4_1m_Test24_dx3cm
$RUNFDSMPI 16 $INDIR Sandia_CH4_1m_Test24_dx1p5cm
$RUNFDSMPI 16 $INDIR Sandia_H2_1m_Test35_dx6cm
$RUNFDSMPI 16 $INDIR Sandia_H2_1m_Test35_dx3cm
$RUNFDSMPI 16 $INDIR Sandia_H2_1m_Test35_dx1p5cm

echo FDS cases submitted
