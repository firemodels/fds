#!/bin/csh -f
setenv SVNROOT ~/FDS-SMV
setenv FDS $SVNROOT/FDS_Compilation/mpi_intel_linux_32/fds5_mpi_intel_linux_32
set WDIR=$SVNROOT/Validation/VTT
cd $WDIR/Current_Results
mpirun n3 n3 n3 n3 n3 $FDS VTT_01_v5.fds > & VTT_01_v5.err &
mpirun n4 n4 n4 n4 n4 $FDS VTT_02_v5.fds > & VTT_02_v5.err &
mpirun n5 n5 n5 n5 n5 $FDS VTT_03_v5.fds > & VTT_03_v5.err &
cd $WDIR



