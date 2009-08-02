#!/bin/csh -f
setenv SVNROOT ~/FDS-SMV
setenv FDS $SVNROOT/FDS_Compilation/mpi_intel_linux_32/fds5_mpi_intel_linux_32
set WDIR=$SVNROOT/Validation/ATF_Corridors
cd $WDIR/Current_Results
mpirun n8  n8  n9  $FDS ATF_Corridors_050_kW.fds >& ATF_Corridors_050_kW.err &
mpirun n9  n10 n10 $FDS ATF_Corridors_100_kW.fds >& ATF_Corridors_100_kW.err &
mpirun n11 n11 n12 $FDS ATF_Corridors_240_kW.fds >& ATF_Corridors_240_kW.err &
mpirun n12 n13 n13 $FDS ATF_Corridors_250_kW.fds >& ATF_Corridors_250_kW.err &
mpirun n14 n14 n15 $FDS ATF_Corridors_500_kW.fds >& ATF_Corridors_500_kW.err &
mpirun n15 n16 n16 $FDS ATF_Corridors_Mix_kW.fds >& ATF_Corridors_Mix_kW.err &
cd $WDIR


