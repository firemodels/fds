# MPI Scaling Tests

The sub-folder `FDS_Input_Files` contains a set of input files used to test the MPI functionality in FDS. There are two sets of tests. The first, `strong_scaling_test_N.fds`, is a simple simulation that runs for 100 time steps. The integer `N` indicates the number of meshes. `N=001` is a single mesh case run with 1 MPI process (and 1 OpenMP thread). `N=008` is this same case, but divided into 8 meshes and run with 8 MPI processes (and still only 1 OpenMP thread). `N=432` is the same case divided into 432 meshes and run with 432 MPI processes. Becaues there is the same number of grid cells in total, the run times decrease with the increasing number of meshes. Ideally, the 432 mesh case would run 432 times faster than the 1 mesh case, but obviously inefficiencies will not allow this.

The cases labelled `weak_scaling_test_N` consists of `N` 50 by 50 by 50 cell meshes forming a line of connected meshes in the x-direction. Each case is run with `N` MPI processes, and ideally each job should take the same amount of wall clock time.

The script `Run_All.sh` can be used to run the cases, and `Process_All.sh` copies the output to the `out` repository.


