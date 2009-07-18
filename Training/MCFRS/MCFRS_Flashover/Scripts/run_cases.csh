#!/bin/csh -f
cd ../FDS_Input_Files
mpirun n3 n3 n3 n3 ~/bin/fds5_mpi_linux MCFRS_Flashover_00leak.fds  >& MCFRS_Flashover_00leak.err &
mpirun n4 n4 n4 n4 ~/bin/fds5_mpi_linux MCFRS_Flashover_00open.fds  >& MCFRS_Flashover_00open.err &
mpirun n5 n5 n5 n5 ~/bin/fds5_mpi_linux MCFRS_Flashover_00case3.fds  >& MCFRS_Flashover_00case3.err &
mpirun n6 n6 n6 n6 ~/bin/fds5_mpi_linux MCFRS_Flashover_00case1.fds  >& MCFRS_Flashover_00case1.err &
#mpirun n7 n7 n7 n7 ~/bin/fds5_mpi_linux MCFRS_Flashover_00leakv.fds  >& MCFRS_Flashover_00leakv.err &
#mpirun n7 n7 n7 n7 ~/bin/fds5_mpi_linux MCFRS_Flashover_00openv.fds  >& MCFRS_Flashover_00openv.err &

