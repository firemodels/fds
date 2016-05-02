@echo off

echo running command: mpiexec -hosts 1 %COMPUTERNAME% 4 test_mpi_write
mpiexec -hosts 1 %COMPUTERNAME% 4 test_mpi_write