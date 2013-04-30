#!/bin/bash

# openmp_inspect.sh
# Kristopher Overholt
# 4/30/2013

# Perform OpenMP thread checking (locate deadlocks and data races)

export SVNROOT=`pwd`/../..
source /opt/intel/inspector_xe/inspxe-vars.sh quiet

cd $SVNROOT/Verification/Timing_Benchmarks
inspxe-cl -collect ti3 -result-dir r000ti3 -- $SVNROOT/FDS_Compilation/openmp_intel_linux_64_db/fds_openmp_intel_linux_64_db bench2.fds
