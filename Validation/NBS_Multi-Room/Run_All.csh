#!/bin/csh -f
setenv SVNROOT ~/FDS-SMV
setenv FDS $SVNROOT/FDS_Compilation/mpi_intel_linux_32/fds5_mpi_intel_linux_32
set WDIR=$SVNROOT/Validation/NBS_Multi-Room
cd $WDIR/Current_Results
mpirun n8 n8     $FDS NBS_100A.fds >& NBS_100A.err &
mpirun n8 n8     $FDS NBS_100O.fds >& NBS_100O.err &
mpirun n8 n8 n8  $FDS NBS_100Z.fds >& NBS_100Z.err &
cd $WDIR


