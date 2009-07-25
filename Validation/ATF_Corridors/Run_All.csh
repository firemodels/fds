#!/bin/csh -f
setenv SVNROOT ~/FDS-SMV
setenv FDS $SVNROOT/FDS_Compilation/mpi_intel_linux_32/fds5_mpi_intel_linux_32
set WDIR=$SVNROOT/Validation/ATF_Corridors
cd $WDIR/Current_Results
mpirun n2 n2 n3     $FDS ATF_Corridors_050_kW.fds >& ATF_Corridors_050_kW.err &
mpirun n3 n4 n4     $FDS ATF_Corridors_100_kW.fds >& ATF_Corridors_100_kW.err &
mpirun n5 n5 n8     $FDS ATF_Corridors_250_kW.fds >& ATF_Corridors_250_kW.err &
mpirun n8 n9 n9     $FDS ATF_Corridors_500_kW.fds >& ATF_Corridors_500_kW.err &
mpirun n10 n10 n12  $FDS ATF_Corridors_Mix_kW.fds >& ATF_Corridors_Mix_kW.err &
cd $WDIR


