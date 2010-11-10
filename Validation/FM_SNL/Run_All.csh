#!/bin/csh -f
setenv SVNROOT ~/FDS-SMV
setenv FDS $SVNROOT/FDS_Compilation/mpi_intel_linux_32/fds5_mpi_intel_linux_32
set WDIR=$SVNROOT/Validation/FM_SNL
cd $WDIR/Current_Results
mpirun n1 n1 n1 n1 n1      $FDS FM_SNL_04.fds > & FM_SNL_04.err &
mpirun n2 n2 n2 n2 n2      $FDS FM_SNL_05.fds > & FM_SNL_05.err &
mpirun n3 n3 n3 n3 n3      $FDS FM_SNL_21.fds > & FM_SNL_21.err &
cd $WDIR

