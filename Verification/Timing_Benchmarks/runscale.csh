#!/bin/csh -f

# specify where to run benchmark

set node=n3

# shouldn't have to edit any lines below

set mpi1=$node
set mpi2=$node $node
set mpi3=$node $node $node $node
set mpi4=$node $node $node $node $node $node $node $node
mpirun $mpi1 ~/bin/fds5_mpi_intel_linux_32 scale1.fds >& scale1.err
mpirun $mpi2 ~/bin/fds5_mpi_intel_linux_32 scale2.fds >& scale2.err
mpirun $mpi3 ~/bin/fds5_mpi_intel_linux_32 scale3.fds >& scale3.err
mpirun $mpi4 ~/bin/fds5_mpi_intel_linux_32 scale4.fds >& scale4.err
