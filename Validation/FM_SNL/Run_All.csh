#!/bin/csh -f
setenv SVNROOT ~/FDS-SMV
setenv FDS $SVNROOT/FDS_Compilation/mpi_intel_linux_32/fds5_mpi_intel_linux_32
set WDIR=$SVNROOT/Validation/FM_SNL
cd $WDIR/Current_Results
mpirun n8 n8 n9 n9 n9      $FDS FM_SNL_04_v5.fds > & FM_SNL_04_v5.err &
mpirun n10 n10 n11 n11 n11 $FDS FM_SNL_05_v5.fds > & FM_SNL_05_v5.err &
mpirun n12 n12 n13 n13 n13 $FDS FM_SNL_21_v5.fds > & FM_SNL_21_v5.err &
cd $WDIR

