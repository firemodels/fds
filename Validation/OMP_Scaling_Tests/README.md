# OpenMP Scaling Tests

The sub-folders `FDS_Input_Files` contains a set of input files used to test the OpenMP functionality in FDS. There are two sets of tests. The first, `openmp_test64n.fds`, is a single mesh case with 64 by 64 by 64 grid cells. It runs for 10~s of simulation time. The case `openmp_test64a.fds` is run with one OpenMP thread, `openmp_test64b.fds` is run with two, and so on up to `openmp_test64h.fds` which is run with eight. The cases are all run on a single dedicated node with no other jobs running. 

The cases labelled `test128n.fds` are similar, but with a mesh of 128 by 128 by 128 grid cells. 

The script `Run_All.sh` can be used to run the cases, and `Process_All.sh` copies the output to the `out` repository.


