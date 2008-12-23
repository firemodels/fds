#!/bin/csh -f
cd ../FDS_Input_Files
mpirun n6 n6 n6 n6 ~/bin/fds5_mpi_linux MCFRS_Flashover_00leak.fds  >& MCFRS_Flashover_00leak.err &
mpirun n6 n6 n6 n6 ~/bin/fds5_mpi_linux MCFRS_Flashover_00open.fds  >& MCFRS_Flashover_00open.err &

