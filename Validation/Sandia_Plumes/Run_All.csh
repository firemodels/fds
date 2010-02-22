#!/bin/csh -f
setenv SVNROOT ~/FDS-SMV
setenv FDS $SVNROOT/FDS_Compilation/mpi_intel_linux_32/fds5_mpi_intel_linux_32
set WDIR=$SVNROOT/Validation/Sandia_Plumes
cd $WDIR/Current_Results
mpirun -np 16 $FDS Sandia_He_1m_dx3cm.fds >& Sandia_He_1m_dx3cm.err &
mpirun -np 16 $FDS Sandia_He_1m_dx1p5cm.fds >& Sandia_He_1m_dx1p5cm.err &
cd $WDIR
echo FDS cases submitted
