#!/bin/csh -f
setenv SVNROOT ~/FDS-SMV
setenv FDS $SVNROOT/FDS_Compilation/mpi_intel_linux_32/fds5_mpi_intel_linux_32
set WDIR=$SVNROOT/Validation/ATF_Corridors
cd $WDIR/Current_Results
mpirun n17 n17 n17 $FDS ATF_Corridors_050_kW.fds >& ATF_Corridors_050_kW.err &
mpirun n18 n18 n18 $FDS ATF_Corridors_100_kW.fds >& ATF_Corridors_100_kW.err &
mpirun n19 n19 n19 $FDS ATF_Corridors_240_kW.fds >& ATF_Corridors_240_kW.err &
mpirun n20 n20 n20 $FDS ATF_Corridors_250_kW.fds >& ATF_Corridors_250_kW.err &
mpirun n21 n21 n21 $FDS ATF_Corridors_500_kW.fds >& ATF_Corridors_500_kW.err &
mpirun n22 n22 n22 $FDS ATF_Corridors_Mix_kW.fds >& ATF_Corridors_Mix_kW.err &
cd $WDIR


