# Compiling FDS with GNU Fortran, OpenMPI and Intel Performance Libraries in Ubuntu Linux

This tutorial should help the advanced users who are interested in compiling the FDS source code with all capabilities, using the latest GNU Fortran (7.2) distribution, OpenMPI and linking to Intel Performance libraries. This will provide a locally built FDS with all capabilites found in the bundle.

FDS GitHub Wiki
- https://github.com/firemodels/fds/wiki/Compiling-FDS-with-GNU-Fortran,-OpenMPI-and-Intel-Performance-Libraries-in-Ubuntu-Linux

## Run Ubuntu 14.04.5 LTS

> vagrant up

> lsb_release -a
No LSB modules are available.
Distributor ID:	Ubuntu
Description:	Ubuntu 14.04.5 LTS
Release:	14.04
Codename:	trusty


## Install Build-Essential

> sudo apt-get install build-essential

gcc --version
gcc (Ubuntu 4.8.5-2ubuntu1~14.04.1) 4.8.5




## Install GNU Fortran

vagrant@vagrant-ubuntu-trusty-64:~$ sudo add-apt-repository ppa:ubuntu-toolchain-r/test
 Toolchain test builds; see https://wiki.ubuntu.com/ToolChain

 More info: https://launchpad.net/~ubuntu-toolchain-r/+archive/ubuntu/test
Press [ENTER] to continue or ctrl-c to cancel adding it

gpg: keyring `/tmp/tmpgve5rsn_/secring.gpg' created
gpg: keyring `/tmp/tmpgve5rsn_/pubring.gpg' created
gpg: requesting key BA9EF27F from hkp server keyserver.ubuntu.com
gpg: /tmp/tmpgve5rsn_/trustdb.gpg: trustdb created
gpg: key BA9EF27F: public key "Launchpad Toolchain builds" imported
gpg: Total number processed: 1
gpg:               imported: 1  (RSA: 1)
OK


vagrant@vagrant-ubuntu-trusty-64:~$ sudo apt install gfortran-7
Reading package lists... Done
Building dependency tree       
Reading state information... Done
The following extra packages will be installed:
  cpp-7 gcc-7 gcc-7-base libasan4 libatomic1 libcc1-0 libcilkrts5 libgcc-7-dev
  libgcc1 libgfortran-7-dev libgfortran4 libgomp1 libisl15 libitm1 liblsan0
  libmpfr4 libmpx2 libquadmath0 libstdc++6 libtsan0 libubsan0
Suggested packages:
  gcc-7-locales gcc-7-multilib gcc-7-doc libgcc1-dbg libgomp1-dbg libitm1-dbg
  libatomic1-dbg libasan4-dbg liblsan0-dbg libtsan0-dbg libubsan0-dbg
  libcilkrts5-dbg libmpx2-dbg libquadmath0-dbg gfortran-7-multilib
  gfortran-7-doc libgfortran4-dbg libcoarrays-dev
The following NEW packages will be installed:
  cpp-7 gcc-7 gcc-7-base gfortran-7 libasan4 libcc1-0 libcilkrts5 libgcc-7-dev
  libgfortran-7-dev libgfortran4 libisl15 liblsan0 libmpx2 libubsan0
The following packages will be upgraded:
  libatomic1 libgcc1 libgomp1 libitm1 libmpfr4 libquadmath0 libstdc++6
  libtsan0
8 upgraded, 14 newly installed, 0 to remove and 8 not upgraded.
Need to get 105 MB of archives.
After this operation, 492 MB of additional disk space will be used.
Do you want to continue? [Y/n] y
...
...
...
Setting up cpp-7 (7.2.0-1ubuntu1~14.04) ...
Setting up libgcc-7-dev:amd64 (7.2.0-1ubuntu1~14.04) ...
Setting up gcc-7 (7.2.0-1ubuntu1~14.04) ...
Setting up libgfortran-7-dev:amd64 (7.2.0-1ubuntu1~14.04) ...
Setting up gfortran-7 (7.2.0-1ubuntu1~14.04) ...
Processing triggers for libc-bin (2.19-0ubuntu6.13) ...



## Compile OpenMPI 3.X

https://www.open-mpi.org/faq/?category=building#easy-build

> gunzip -c openmpi-3.0.0.tar.gz | tar xf -
> cd openmpi-3.0.0/
> ./configure FC=gfortran-7 CC=gcc-7 --prefix=/shared/openmpi_64 --enable-mpirun-prefix-by-default --enable-mpi-fortran --enable-static --disable-shared


> make -j 2

Reduce build time `make -j 4` 
> man make

       -j [jobs], --jobs[=jobs]
            Specifies  the  number  of  jobs  (commands) to run simultaneously.  If there is more than one -j option, the last one is effective.  If the -j option is given without an
            argument, make will not limit the number of jobs that can run simultaneously.

	Some versions of make support parallel builds. The example above shows GNU make's "-j" option, which specifies how many compile processes may be executing any any given time. We, the Open MPI Team, have found that doubling or quadrupling the number of processors in a machine can significantly speed up an Open MPI compile (since compiles tend to be much more IO bound than CPU bound).

FAQ - Building Open MPI
- https://www.open-mpi.org/faq/?category=building#vpath-parallel-build

> sudo make install

...

Libraries have been installed in:
   /shared/openmpi_64/lib

If you ever happen to want to link against installed libraries
in a given directory, LIBDIR, you must either use libtool, and
specify the full pathname of the library, or use the '-LLIBDIR'
flag during linking and do at least one of the following:
   - add LIBDIR to the 'LD_LIBRARY_PATH' environment variable
     during execution
   - add LIBDIR to the 'LD_RUN_PATH' environment variable
     during linking
   - use the '-Wl,-rpath -Wl,LIBDIR' linker flag
   - have your system administrator add LIBDIR to '/etc/ld.so.conf'

See any operating system documentation about shared libraries for
more information, such as the ld(1) and ld.so(8) manual pages.

vagrant@vagrant-ubuntu-trusty-64:/vagrant/Build/mpi_gnu_linux_64$ which mpirun
/shared/openmpi_64/bin/mpirun
vagrant@vagrant-ubuntu-trusty-64:/vagrant/Build/mpi_gnu_linux_64$ mpirun
--------------------------------------------------------------------------
mpirun could not find anything to do.

It is possible that you forgot to specify how many processes to run
via the "-np" argument.
--------------------------------------------------------------------------
vagrant@vagrant-ubuntu-trusty-64:/vagrant/Build/mpi_gnu_linux_64$ mpirun --version
mpirun (Open MPI) 3.0.0


## Compile FDS sources

Make these directories available to your environment:

export MPIDIST=/shared/openmpi_64
export PATH=$MPIDIST/bin:$PATH
export LD_LIBRARY_PATH=$MPIDIST/lib:$LD_LIBRARY_PATH




vagrant@vagrant-ubuntu-trusty-64:/vagrant/Build/mpi_gnu_linux_64$ ./fds_mpi_gnu_linux_64 

 Fire Dynamics Simulator

 Current Date     : November 22, 2017  12:31:16
 Version          : FDS 6.6.0
 Revision         : -
 Revision Date    : 
 Compiler         : unknown
 Compilation Date : Nov 22, 2017  12:05:17

 MPI Enabled;    Number of MPI Processes:       1
 OpenMP Enabled; Number of OpenMP Threads:      2

 MPI version: 3.1
 MPI library version: Open MPI v3.0.0, package: Open MPI vagrant@vagrant-ubuntu-trusty-64 Distribution, ident: 3.0.0, repo rev: v3.0.0, Sep 12, 2017

 Consult FDS Users Guide Chapter, Running FDS, for further instructions.

 Hit Enter to Escape...




## Compilers

GNU Fortran 8.0
- https://www.scivision.co/install-latest-gfortran-on-ubuntu/


## CI as a service

Build FDS on Linux
- Travis CI

Travis CI - Install dependencies
 - https://docs.travis-ci.com/user/installing-dependencies/#Installing-Packages-from-a-custom-APT-repository



## Useful tools 

Michael Hirsch
- https://www.scivision.co/install-linuxbrew/