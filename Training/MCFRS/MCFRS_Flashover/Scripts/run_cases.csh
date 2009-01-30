#!/bin/csh -f
cd ../FDS_Input_Files
mpirun n0 n0 n0 n0 ~/bin/fds5_mpi_linux MCFRS_Flashover_00leak.fds  >& MCFRS_Flashover_00leak.err &
mpirun n1 n1 n1 n1 ~/bin/fds5_mpi_linux MCFRS_Flashover_00open.fds  >& MCFRS_Flashover_00open.err &
mpirun n2 n2 n2 n2 ~/bin/fds5_mpi_linux MCFRS_Flashover_00case3.fds  >& MCFRS_Flashover_00case3.err &
#mpirun n7 n7 n7 n7 ~/bin/fds5_mpi_linux MCFRS_Flashover_00leakv.fds  >& MCFRS_Flashover_00leakv.err &
#mpirun n7 n7 n7 n7 ~/bin/fds5_mpi_linux MCFRS_Flashover_00openv.fds  >& MCFRS_Flashover_00openv.err &

