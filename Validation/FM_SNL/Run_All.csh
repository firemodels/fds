#!/bin/csh -f
setenv SVNROOT ~/FDS-SMV
setenv FDS $SVNROOT/FDS_Compilation/mpi_intel_linux_32/fds5_mpi_intel_linux_32
set WDIR=$SVNROOT/Validation/FM_SNL
cd $WDIR/Current_Results
mpirun n0 n0 n0 n0 n0      $FDS FM_SNL_04_v5.fds > & FM_SNL_04_v5.err &
mpirun n1 n1 n1 n1 n1      $FDS FM_SNL_05_v5.fds > & FM_SNL_05_v5.err &
mpirun n2 n2 n2 n2 n2      $FDS FM_SNL_21_v5.fds > & FM_SNL_21_v5.err &
cd $WDIR

