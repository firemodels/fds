#!/bin/csh -f
setenv SVNROOT ~/FDS-SMV
setenv FDS $SVNROOT/FDS_Compilation/mpi_intel_linux_32/fds5_mpi_intel_linux_32
set WDIR=$SVNROOT/Validation/UL_SFPE
cd $WDIR/Current_Results
mpirun n11 n11     $FDS UL_SFPE_10.fds >& UL_SFPE_10.err &
mpirun n13 n13     $FDS UL_SFPE_15.fds >& UL_SFPE_15.err &
mpirun n14 n14     $FDS UL_SFPE_20.fds >& UL_SFPE_20.err &
mpirun n15 n15     $FDS UL_SFPE_25.fds >& UL_SFPE_25.err &
mpirun n16 n16     $FDS UL_SFPE_35.fds >& UL_SFPE_35.err &
mpirun n17 n17     $FDS UL_SFPE_40.fds >& UL_SFPE_40.err &
cd $WDIR


