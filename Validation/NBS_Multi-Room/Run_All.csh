#!/bin/csh -f
setenv SVNROOT ~/FDS-SMV
setenv FDS $SVNROOT/FDS_Compilation/mpi_intel_linux_32/fds5_mpi_intel_linux_32
set WDIR=$SVNROOT/Validation/NBS_Multi-Room
cd $WDIR/Current_Results
mpirun n23 n23     $FDS NBS_100A_v5.fds >& NBS_100A_v5.err &
mpirun n14 n14     $FDS NBS_100O_v5.fds >& NBS_100O_v5.err &
mpirun n15 n15 n16 $FDS NBS_100Z_v5.fds >& NBS_100Z_v5.err &
cd $WDIR


