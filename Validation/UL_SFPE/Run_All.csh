#!/bin/csh -f
setenv SVNROOT ~/FDS-SMV
setenv FDS $SVNROOT/FDS_Compilation/mpi_intel_linux_32/fds5_mpi_intel_linux_32
set WDIR=$SVNROOT/Validation/UL_SFPE
cd $WDIR/Current_Results
mpirun n15 n15     $FDS UL_SFPE_10.fds >& UL_SFPE_10.err &
mpirun n11 n11     $FDS UL_SFPE_15.fds >& UL_SFPE_15.err &
mpirun n10 n10     $FDS UL_SFPE_20.fds >& UL_SFPE_20.err &
mpirun n9  n9      $FDS UL_SFPE_25.fds >& UL_SFPE_25.err &
mpirun n4  n7      $FDS UL_SFPE_35.fds >& UL_SFPE_35.err &
mpirun n6  n6      $FDS UL_SFPE_40.fds >& UL_SFPE_40.err &
cd $WDIR


