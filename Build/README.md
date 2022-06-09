## Compiling FDS

The folders above are used for compiling FDS with different versions of MPI (impi is for Intel MPI; ompi is for Open MPI), different Fortran compilers (intel, gnu), different operating systems (Windows, Linux, Mac OS X), and mode (db means debug, dv means low level of optimization, no suffix indicates release version). The attribute `openmp` indicates that the executable has been built with OpenMP directives. For the Intel Fortran compiler, it was found that compiling with OpenMP incurs an overhead tax of 20%. All builds that use the Gnu Fortran compiler include OpenMP.

All folders contain a single bash script called `make_fds.sh` that invokes the same `makefile`. All the build targets are listed in the `makefile`. 

The links to the Intel MKL libraries were created using a [tool provided by Intel](https://www.intel.com/content/www/us/en/developer/tools/oneapi/onemkl-link-line-advisor.html#gs.xoagc3).

If you find a particular set of options that is no longer appropriate, let us know, preferably via a pull request for the makefile.



