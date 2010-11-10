#!/bin/csh -f
setenv SVNROOT ~/FDS-SMV
setenv FDS $SVNROOT/FDS_Compilation/mpi_intel_linux_32/fds5_mpi_intel_linux_32
set WDIR=$SVNROOT/Validation/VTT
cd $WDIR/Current_Results
mpirun n6 n6 n6 n6 n6 $FDS VTT_01.fds > & VTT_01.err &
mpirun n4 n4 n4 n4 n4 $FDS VTT_02.fds > & VTT_02.err &
mpirun n5 n5 n5 n5 n5 $FDS VTT_03.fds > & VTT_03.err &
cd $WDIR



