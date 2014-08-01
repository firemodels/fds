#!/bin/bash

# inspect_openmp.sh
# Kristopher Overholt
# 4/30/2013

# Perform OpenMP thread checking (locate deadlocks and data races)

export SVNROOT=`pwd`/../..
source /opt/intel/inspector_xe/inspxe-vars.sh quiet

TARGET=$SVNROOT/FDS_Compilation/intel_linux_64_inspect

if [ -d inspect_openmp_ti3 ]
then
    rm -r inspect_openmp_ti3
fi

cd $TARGET
make -f ../makefile clean
./make_fds.sh

cd $SVNROOT/Verification/Timing_Benchmarks
export OMP_NUM_THREADS=2
inspxe-cl -collect ti3 -knob scope=normal \
          -result-dir $SVNROOT/Utilities/Scripts/inspect_openmp_ti3 \
          -search-dir src=$SVNROOT/FDS_Source \
          -- $SVNROOT/FDS_Compilation/intel_linux_64_inspect/fds_intel_linux_64_inspect bench2.fds
