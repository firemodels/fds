#!/bin/csh -f
setenv SVNROOT ~/FDS-SMV
setenv FDS $SVNROOT/FDS_Compilation/mpi_intel_linux_32/fds5_mpi_intel_linux_32
set WDIR=$SVNROOT/Validation/VTT
cd $WDIR/Current_Results
mpirun n17 n17 n17 n18 n18 $FDS VTT_01_v5.fds > & VTT_01_v5.err &
mpirun n19 n19 n20 n20 n20 $FDS VTT_02_v5.fds > & VTT_02_v5.err &
mpirun n21 n21 n22 n22 n22 $FDS VTT_03_v5.fds > & VTT_03_v5.err &
cd $WDIR



