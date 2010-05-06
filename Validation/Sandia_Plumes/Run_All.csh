#!/bin/csh -f
setenv SVNROOT ~/FDS-SMV
setenv FDS $SVNROOT/FDS_Compilation/mpi_intel_linux_32/fds5_mpi_intel_linux_32
set WDIR=$SVNROOT/Validation/Sandia_Plumes
cd $WDIR/Current_Results
mpirun -np 16 $FDS Sandia_He_1m_dx6cm.fds   >& Sandia_He_1m_dx6cm.err &
mpirun -np 16 $FDS Sandia_He_1m_dx3cm.fds   >& Sandia_He_1m_dx3cm.err &
mpirun -np 16 $FDS Sandia_He_1m_dx1p5cm.fds >& Sandia_He_1m_dx1p5cm.err &
mpirun -np 16 $FDS Sandia_CH4_1m_Test14_dx6cm.fds   >& Sandia_CH4_1m_Test14_dx6cm.err &
mpirun -np 16 $FDS Sandia_CH4_1m_Test14_dx3cm.fds   >& Sandia_CH4_1m_Test14_dx3cm.err &
mpirun -np 16 $FDS Sandia_CH4_1m_Test14_dx1p5cm.fds >& Sandia_CH4_1m_Test14_dx1p5cm.err &
mpirun -np 16 $FDS Sandia_CH4_1m_Test17_dx6cm.fds   >& Sandia_CH4_1m_Test17_dx6cm.err &
mpirun -np 16 $FDS Sandia_CH4_1m_Test17_dx3cm.fds   >& Sandia_CH4_1m_Test17_dx3cm.err &
mpirun -np 16 $FDS Sandia_CH4_1m_Test17_dx1p5cm.fds >& Sandia_CH4_1m_Test17_dx1p5cm.err &
mpirun -np 16 $FDS Sandia_CH4_1m_Test24_dx6cm.fds   >& Sandia_CH4_1m_Test24_dx6cm.err &
mpirun -np 16 $FDS Sandia_CH4_1m_Test24_dx3cm.fds   >& Sandia_CH4_1m_Test24_dx3cm.err &
mpirun -np 16 $FDS Sandia_CH4_1m_Test24_dx1p5cm.fds >& Sandia_CH4_1m_Test24_dx1p5cm.err &
mpirun -np 16 $FDS Sandia_H2_1m_Test35_dx6cm.fds   >& Sandia_H2_1m_Test35_dx6cm.err &
mpirun -np 16 $FDS Sandia_H2_1m_Test35_dx3cm.fds   >& Sandia_H2_1m_Test35_dx3cm.err &
mpirun -np 16 $FDS Sandia_H2_1m_Test35_dx1p5cm.fds >& Sandia_H2_1m_Test35_dx1p5cm.err &

cd $WDIR
echo FDS cases submitted
