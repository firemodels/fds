# Compiling FDS

The folders above are used for compiling FDS with different versions of MPI (impi is for Intel MPI; mpi is for Open MPI), different Fortran compilers (Intel, Gnu, PGI), different operating systems (Windows, Linux, Mac OS X), and mode (db means debug, dv means low level of optimization, no suffix indicates release version). There are also some special builds used for analysis such as `vtune`, `inspect`, `advise`, and `trace`. Details on these are found in the notes under the respective directories.

All folders contain a single bash script called `make_fds.sh` that invokes the same `makefile`. All the build targets are listed in the `makefile`. 

If you find a particular set of options that is no longer appropriate, let us know, preferably via a pull request for the makefile.
