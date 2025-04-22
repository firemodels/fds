@echo off
set I_MPI_ROOT=%~dp0\mpi
set PATH=%I_MPI_ROOT%;%PATH%
set IN_CMDFDS=1
set MPIEXEC_PORT_RANGE=
set MPICH_PORT_RANGE=

title FDS
echo.
echo type helpfds for help on running fds

