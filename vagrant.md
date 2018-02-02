# Compiling FDS with GNU Fortran, OpenMPI and Intel Performance Libraries in Ubuntu Linux

This tutorial should help the advanced users who are interested in compiling the FDS source code with all capabilities, using the latest GNU Fortran (7.2) distribution, OpenMPI and linking to Intel Performance libraries. This will provide a locally built FDS with all capabilites found in the bundle.

FDS GitHub Wiki
- [Compiling FDS with GNU Fortran, OpenMPI and Intel Performance Libraries in Ubuntu Linux](https://github.com/firemodels/fds/wiki/Compiling-FDS-with-GNU-Fortran,-OpenMPI-and-Intel-Performance-Libraries-in-Ubuntu-Linux)

## Spin vagrant box running Ubuntu 14.04.5 LTS (Trusty Tar)

    $ vagrant up

    $ lsb_release -a
    No LSB modules are available.
    Distributor ID:	Ubuntu
    Description:	Ubuntu 14.04.5 LTS
    Release:	14.04
    Codename:	trusty


## Install package build-essential

    $ sudo apt-get install build-essential

    $ gcc --version
    gcc (Ubuntu 4.8.5-2ubuntu1~14.04.1) 4.8.5


## Install GNU Fortran

Add Ubuntu PPA repository

    $ sudo add-apt-repository ppa:ubuntu-toolchain-r/test
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

Install GNU Fortran 7

    $ sudo apt install gfortran-7
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


Verify Fortran installation

    $ gfortran-7 --version
    GNU Fortran (Ubuntu 7.2.0-1ubuntu1~14.04) 7.2.0
    Copyright (C) 2017 Free Software Foundation, Inc.
    This is free software; see the source for copying conditions.  There is NO
    warranty; not even for MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.


## Compile OpenMPI 3.X

Download from official website
    wget https://www.open-mpi.org/software/ompi/v3.0/downloads/openmpi-3.0.0.tar.gz




    $ gunzip -c openmpi-3.0.0.tar.gz | tar xf -
    $ cd openmpi-3.0.0/
    $ ./configure FC=gfortran-7 CC=gcc-7 --prefix=/shared/openmpi_64 --enable-mpirun-prefix-by-default --enable-mpi-fortran --enable-static --disable-shared
    $ make -j 2

Explain installation options

    --prefix=<directory>
      Install Open MPI into the base directory named <directory>.  Hence,
      Open MPI will place its executables in <directory>/bin, its header
      files in <directory>/include, its libraries in <directory>/lib, etc.



    --disable-shared
      By default, Open MPI and OpenSHMEM build shared libraries, and all
      components are built as dynamic shared objects (DSOs). This switch
      disables this default; it is really only useful when used with
      --enable-static.  Specifically, this option does *not* imply
      --enable-static; enabling static libraries and disabling shared
      libraries are two independent options.

    --enable-static
      Build MPI and OpenSHMEM as static libraries, and statically link in
      all components.  Note that this option does *not* imply
      --disable-shared; enabling static libraries and disabling shared
      libraries are two independent options.

    --enable-mpi-fortran(=value)
      By default, Open MPI will attempt to build all 3 Fortran bindings:
      mpif.h, the "mpi" module, and the "mpi_f08" module.  The following
      values are permitted:

        all:        Synonym for "yes".
        yes:        Attempt to build all 3 Fortran bindings; skip
                    any binding that cannot be built (same as
                    --enable-mpi-fortran).
        mpifh:      Build mpif.h support.
        usempi:     Build mpif.h and "mpi" module support.
        usempif08:  Build mpif.h, "mpi" module, and "mpi_f08"
                    module support.
        none:       Synonym for "no".
        no:         Do not build any MPI Fortran support (same as
                    --disable-mpi-fortran).  This is mutually exclusive
                    with building the OpenSHMEM Fortran interface.

RUN-TIME SYSTEM SUPPORT

Note that in both models, invoking mpirun via an absolute path name is equivalent to specifying the --prefix option with a <dir> value equivalent to the directory where mpirun resides, minus its last subdirectory. For example:

    % /usr/local/bin/mpirun ...

is equivalent to

    % mpirun --prefix /usr/local

    --enable-mpirun-prefix-by-default
      This option forces the "mpirun" command to always behave as if
      "--prefix $prefix" was present on the command line (where $prefix is
      the value given to the --prefix option to configure).  This prevents
      most rsh/ssh-based users from needing to modify their shell startup
      files to set the PATH and/or LD_LIBRARY_PATH for Open MPI on remote
      nodes.  Note, however, that such users may still desire to set PATH
      -- perhaps even in their shell startup files -- so that executables
      such as mpicc and mpirun can be found without needing to type long
      path names.  --enable-orterun-prefix-by-default is a synonym for
      this option.

    --noprefix  
      Disable the automatic --prefix behavior 




Reduce build time

    $ make -j 4

Why this works? Explain `make` argument `j`

    $ man make

         -j [jobs], --jobs[=jobs]
              Specifies  the  number  of  jobs  (commands) to run simultaneously.  If there is more than one -j option, the last one is effective.  If the -j option is given without an
              argument, make will not limit the number of jobs that can run simultaneously.

Open MPI teams comment on this

>Some versions of make support parallel builds. The example above shows GNU make's "-j" option, which specifies how many compile processes may be executing any any given time. We, the Open MPI Team, have found that doubling or quadrupling the number of processors in a machine can significantly speed up an Open MPI compile (since compiles tend to be much more IO bound than CPU bound).



    sudo make install

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

Links

- [Building Open MPI](https://www.open-mpi.org/faq/?category=building#easy-build)
- [FAQ - Building Open MPI](https://www.open-mpi.org/faq/?category=building#vpath-parallel-build)

## Compile FDS sources

Make these directories available to your environment:

    export MPIDIST=/shared/openmpi_64
    export PATH=$MPIDIST/bin:$PATH
    export LD_LIBRARY_PATH=$MPIDIST/lib:$LD_LIBRARY_PATH


Go to your updated repo for FDS, you should be able to build the mpi_gnu_linux_64

    $ cd Build/mpi_gnu_linux_64
    $ ./make_fds.sh

Run compiled FDS exectuable

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


## Run FDS as MPI application

    $ /shared/openmpi_64/bin/mpirun -quiet /vagrant/Build/mpi_gnu_linux_64/fds_mpi_gnu_linux_64

     Fire Dynamics Simulator

     Current Date     : November 22, 2017  13:44:27
     Version          : FDS 6.6.0
     Revision         : -
     Revision Date    :
     Compiler         : unknown
     Compilation Date : Nov 22, 2017  12:05:17

     MPI Enabled;    Number of MPI Processes:       2
     OpenMP Enabled; Number of OpenMP Threads:      1

     MPI version: 3.1
     MPI library version: Open MPI v3.0.0, package: Open MPI vagrant@vagrant-ubuntu-trusty-64 Distribution, ident: 3.0.0, repo rev: v3.0.0, Sep 12, 2017

     Consult FDS Users Guide Chapter, Running FDS, for further instructions.

     Hit Enter to Escape...


## Explore required dependencies of `FDS` executable

Use `ldd` tool to print shared libraries dependencies.


    $ ldd fds_mpi_gnu_linux_64 
        linux-vdso.so.1 =>  (0x00007ffc427ec000)
        libgfortran.so.4 => /usr/lib/x86_64-linux-gnu/libgfortran.so.4 (0x00007f8ab015d000)
        libm.so.6 => /lib/x86_64-linux-gnu/libm.so.6 (0x00007f8aafe57000)
        libdl.so.2 => /lib/x86_64-linux-gnu/libdl.so.2 (0x00007f8aafc53000)
        librt.so.1 => /lib/x86_64-linux-gnu/librt.so.1 (0x00007f8aafa4b000)
        libutil.so.1 => /lib/x86_64-linux-gnu/libutil.so.1 (0x00007f8aaf848000)
        libgomp.so.1 => /usr/lib/x86_64-linux-gnu/libgomp.so.1 (0x00007f8aaf619000)
        libgcc_s.so.1 => /lib/x86_64-linux-gnu/libgcc_s.so.1 (0x00007f8aaf402000)
        libpthread.so.0 => /lib/x86_64-linux-gnu/libpthread.so.0 (0x00007f8aaf1e4000)
        libc.so.6 => /lib/x86_64-linux-gnu/libc.so.6 (0x00007f8aaee1b000)
        libquadmath.so.0 => /usr/lib/x86_64-linux-gnu/libquadmath.so.0 (0x00007f8aaebdc000)
        /lib64/ld-linux-x86-64.so.2 (0x00007f8ab0530000)

Given the dependency of `libgfortran.so.4` a user might install `libgfortran`.

    gottfried@gottfried-SVE14A2X1EH:~/repos/fds/Build/mpi_gnu_linux_64$ sudo apt-get install gfortran
    ...   
    The following additional packages will be installed:
      gfortran-7 libgfortran-7-dev libgfortran4
    Suggested packages:
      gfortran-multilib gfortran-doc gfortran-7-multilib gfortran-7-doc
      libgfortran4-dbg libcoarrays-dev
    The following NEW packages will be installed:
      gfortran gfortran-7 libgfortran-7-dev libgfortran4
    0 upgraded, 4 newly installed, 0 to remove and 14 not upgraded.
    Need to get 8.115 kB of archives.
    After this operation, 30,0 MB of additional disk space will be used.
    Do you want to continue? [Y/n] y
    Get:1 http://de.archive.ubuntu.com/ubuntu artful/main amd64 libgfortran4 amd64 7.2.0-8ubuntu3 [493 kB]
    Get:2 http://de.archive.ubuntu.com/ubuntu artful/main amd64 libgfortran-7-dev amd64 7.2.0-8ubuntu3 [530 kB]
    Get:3 http://de.archive.ubuntu.com/ubuntu artful/main amd64 gfortran-7 amd64 7.2.0-8ubuntu3 [7.090 kB]
    Get:4 http://de.archive.ubuntu.com/ubuntu artful/main amd64 gfortran amd64 4:7.2.0-1ubuntu1 [1.278 B]
    Fetched 8.115 kB in 4s (1.679 kB/s)  
    Selecting previously unselected package libgfortran4:amd64.
    (Reading database ... 409697 files and directories currently installed.)
    Preparing to unpack .../libgfortran4_7.2.0-8ubuntu3_amd64.deb ...
    Unpacking libgfortran4:amd64 (7.2.0-8ubuntu3) ...
    Selecting previously unselected package libgfortran-7-dev:amd64.
    Preparing to unpack .../libgfortran-7-dev_7.2.0-8ubuntu3_amd64.deb ...
    Unpacking libgfortran-7-dev:amd64 (7.2.0-8ubuntu3) ...
    Selecting previously unselected package gfortran-7.
    Preparing to unpack .../gfortran-7_7.2.0-8ubuntu3_amd64.deb ...
    Unpacking gfortran-7 (7.2.0-8ubuntu3) ...
    Selecting previously unselected package gfortran.
    Preparing to unpack .../gfortran_4%3a7.2.0-1ubuntu1_amd64.deb ...
    Unpacking gfortran (4:7.2.0-1ubuntu1) ...
    Setting up libgfortran4:amd64 (7.2.0-8ubuntu3) ...
    Processing triggers for libc-bin (2.26-0ubuntu2.1) ...
    Setting up libgfortran-7-dev:amd64 (7.2.0-8ubuntu3) ...
    Processing triggers for man-db (2.7.6.1-2) ...
    Setting up gfortran-7 (7.2.0-8ubuntu3) ...
    Setting up gfortran (4:7.2.0-1ubuntu1) ...
    update-alternatives: using /usr/bin/gfortran to provide /usr/bin/f95 (f95) in auto mode
    update-alternatives: using /usr/bin/gfortran to provide /usr/bin/f77 (f77) in auto mode


## Deploy FDS/MPI to clusters

- Disable the automatic prefix behaviour to provide OpenMPI in a custom directory

  $ /absolute/path/to/mpirun --noprefix mpi_app

References
- [mpirun - OpenMPI docs](https://www.open-mpi.org/doc/v3.0/man1/mpirun.1.php) 





## More on compilers

GNU Fortran 8.0
- https://www.scivision.co/install-latest-gfortran-on-ubuntu/


## Continuous integration as a service

Build FDS on Linux
- Travis CI - free for open source projects (linux & MacOSX)
- AppVeyor - free for open source projects (Windows)
- Visual Studio Team Services - (Linux, MacOSX, Windows)

Travis CI

 - [Install dependencies](https://docs.travis-ci.com/user/installing-dependencies/#Installing-Packages-from-a-custom-APT-repository)
 - [Embed badges](https://docs.travis-ci.com/user/status-images/)

## Useful tooling

Michael Hirsch
- https://www.scivision.co/install-linuxbrew/

Fortran building system
- https://github.com/szaghi/FoBiS
