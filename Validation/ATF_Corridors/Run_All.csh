#!/bin/csh -f
setenv SVNROOT ~/FDS-SMV
setenv FDS $SVNROOT/FDS_Compilation/mpi_intel_linux_32/fds5_mpi_intel_linux_32
set WDIR=$SVNROOT/Validation/ATF_Corridors
cd $WDIR/Current_Results
mpirun n27 n27 n27 $FDS ATF_Corridors_050_kW.fds >& ATF_Corridors_050_kW.err &
mpirun n28 n28 n28 $FDS ATF_Corridors_100_kW.fds >& ATF_Corridors_100_kW.err &
mpirun n29 n29 n29 $FDS ATF_Corridors_240_kW.fds >& ATF_Corridors_240_kW.err &
mpirun n30 n30 n30 $FDS ATF_Corridors_250_kW.fds >& ATF_Corridors_250_kW.err &
mpirun n31 n31 n31 $FDS ATF_Corridors_500_kW.fds >& ATF_Corridors_500_kW.err &
mpirun n32 n32 n32 $FDS ATF_Corridors_Mix_kW.fds >& ATF_Corridors_Mix_kW.err &
cd $WDIR


