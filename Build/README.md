## Compiling FDS

The folders above are used for compiling FDS with different versions of MPI (impi is for Intel MPI; mpi is for Open MPI), different Fortran compilers (Intel, Gnu, PGI), different operating systems (Windows, Linux, Mac OS X), and mode (db means debug, dv means low level of optimization, no suffix indicates release version). There are also some special builds used for analysis such as `vtune`, `inspect`, `advise`, and `trace`. Details on these are found in the notes under the respective directories.

All folders contain a single bash script called `make_fds.sh` that invokes the same `makefile`. All the build targets are listed in the `makefile`. 

If you find a particular set of options that is no longer appropriate, let us know, preferably via a pull request for the makefile.

## Intel Cluster Checker

The Intel Cluster Checker is a diagnostic tool that checks the overall health of your compute cluster. To use it, first consult the [User's Guide](https://software.intel.com/en-us/cluster-checker-user-guide-2019-beta). In brief, do the following:

   1. Install Intel Cluster Checker and run `source /opt/intel19/clck_latest/clckvars.sh`. Your path might be slightly different.
   2. Create a `nodefile`, which is just a text file with a list of the cluster node names, one per line.
   3. Run `clck -f nodefile`

