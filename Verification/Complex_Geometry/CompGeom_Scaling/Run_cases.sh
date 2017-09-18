#!/bin/bash

# This script runs a set of Validation Cases on a Linux machine with a batch queuing system.
# Executable and base directory must be provided:
export BASEDIR=`pwd`/../../../..
export QFDSDIR=$BASEDIR/fds/Utilities/Scripts
export EXECDIR=$BASEDIR/fds/Build/mpi_intel_linux_64ib
export EXECBIN=/fds_mpi_intel_linux_64ib
export QFDS=$QFDSDIR"/qfds.sh -f "$BASEDIR" -e "$EXECDIR$EXECBIN


# Scaling with increasing number of triangles:
$QFDS -p 1 -o 1 compgeom_scale_64x64_1mesh_320T.fds
$QFDS -p 1 -o 1 compgeom_scale_64x64_1mesh_1280T.fds
$QFDS -p 1 -o 1 compgeom_scale_64x64_1mesh_5120T.fds
$QFDS -p 1 -o 1 compgeom_scale_64x64_1mesh_20480T.fds
$QFDS -p 1 -o 1 compgeom_scale_64x64_1mesh_81920T.fds
$QFDS -p 1 -o 1 compgeom_scale_64x64_1mesh_327680T.fds

# Scaling with increasing number of cells:
$QFDS -p 1 -o 1 compgeom_scale_32x32_1mesh_5120T.fds
$QFDS -p 1 -o 1 compgeom_scale_128x128_1mesh_5120T.fds
$QFDS -p 1 -o 1 compgeom_scale_256x256_1mesh_5120T.fds
