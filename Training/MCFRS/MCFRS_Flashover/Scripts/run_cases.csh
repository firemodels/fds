#!/bin/csh -f
cd ../FDS_Input_Files
mpirun n9 n9 n9 n9 ~/bin/fds5_mpi_linux MCFRS_Flashover_00leak.fds  >& MCFRS_Flashover_00leak.err &
mpirun n9 n9 n9 n9 ~/bin/fds5_mpi_linux MCFRS_Flashover_00open.fds  >& MCFRS_Flashover_00open.err &

