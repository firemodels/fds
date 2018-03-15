# Compiling FDS with GNU Fortran, Open MPI and Intel Performance Libraries in Ubuntu Linux

The basic idea is to create a build environment to compile FDS source code in an automated way - reproducible, version controlled and fast.
Vagrant can provide this to every developer. A file `Vagrantfile` to create and provision a virtual machine is now added to the source code repository in top-level folder.

With regard to running FDS it is a good idea to isolate an end-user environment used for production jobs from an environment capable to compile FDS Fortran source code. End-users just install the official bundle and usually only care about system administration when they have issues. When support end-users reconfigure their environmentto compile nstall software dependencies required


## Why vagrant

If you are a developer, Vagrant will isolate dependencies and their configuration within a single disposable, consistent environment, without sacrificing any of the tools you are used to working with (editors, browsers, debuggers, etc.). Once you or someone else creates a single `Vagrantfile`, you just need to vagrant up and everything is installed and configured for you to work. Other members of your team create their development environments from the same configuration, so whether you are working on Linux, Mac OS X, or Windows, all your team members are running code in the same environment, against the same dependencies, all configured the same way. Say goodbye to "works on my machine" bugs.
The Vagrantfile is meant to be committed to version control with your project, if you use version control. This way, every person working with that project can benefit from Vagrant without any upfront work.

[Vagrant might be for everyone](https://www.vagrantup.com/intro/index.html#for-everyone)

## Install Vagrant 

Download Vagrant
- [Downloads for latest version of Vagrant](https://www.vagrantup.com/downloads.html)
    - Please download the proper binary package for your operating system and architecture. 

Verify installation

```bash
$ vagrant --version
Vagrant 1.9.8
```

>The installer will automatically add vagrant to your system path so that it is available in terminals.

## Install a provider

Vagrant comes with support out of the box for VirtualBox, a free, cross-platform consumer virtualization product.

VirtualBox must be installed on its own prior to using the provider, or the provider will display an error message asking you to install it. VirtualBox can be installed by downloading a package or installer for your operating system and using standard procedures to install that package.

- [Download a package or installer for your operating system](https://www.virtualbox.org/wiki/Downloads)

## Spin up Vagrant box 

The build environment uses Ubuntu 14.04.5 LTS (Trusty Tar) to demonstrate compiling FDS using free and open source software (FOSS). Boot your Vagrant environment. 

    $ vagrant up

## Provisioning

Provisioners in Vagrant allow you to automatically install software, alter configurations, and more on the machine as part of the `vagrant up` process.

Create a build environment
- Install GNU Fortran 7.2
- Build Open MPI 3.0.0
- Install Intel Math Kernel Library
- Compile FDS sources
  - build target `mpi_gnu_linux_64`

## Connect to Vagrant box

 Vagrant runs the virtual machine without a UI. To prove that it is running, you can SSH into the machine

    $ vagrant ssh

## Run FDS executable on Vagrant box

By default, Vagrant shares your project directory (remember, that is the one with the Vagrantfile) to the `/vagrant` directory in your guest machine.

Note that when you `vagrant ssh` into your machine, you're in `/home/vagrant`. `/home/vagrant` is a different directory from the synced `/vagrant` directory.
    
    $ pwd
    /home/vagrant

Run FDS in your home directory

    $ /vagrant/Build/mpi_gnu_linux_64/fds_mpi_gnu_linux_64 case.fds 


## Teardown Vagrant box

How do we clean up our devlopment environment? This will remove all traces of the virtual machine from your system.

    $ vagrant destroy

Or if something goes wrong, just teardown your box and start from scratch with

    $ vagrant destroy
    $ vagrant up

More help on commands to suspend, halt or destroy Vagrant boxes is available on [Vagrant Docs](https://vagrantup.com/intro/getting-started/teardown.html)     

## Support for `GLMAT` pressure solver

To use `GLMAT` pressure solver in your models FDS must be compiled with Intel MKL library. Support for automated setup of Intel MKL in this build environment was added just recently.

Actually there is an issue with the automated way to enable MKL environment. Please use the steps below to set environment for Intel MKL manually.

Connect to Vagrant box

    $ vagrant ssh

Inside Vagrant box

    $ source /opt/intel/mkl/bin/mklvars.sh intel64
    $ source /vagrant/open-mpi-vars.sh
    $ cd /vagrant/Build/mpi_gnu_linux_64
    $ rm -rf *.mod *.o
    $ ./make_fds.sh

>It is possible that you need to delete existing `.mod` and `.o` files in build directory compiled without linking to MKL library.

Verify working Intel MKL environment

> watch out for `DWITH_MKL` flags in compiler and linker commands

```bash
vagrant@vagrant:/vagrant/Build/mpi_gnu_linux_64$ ./make_fds.sh 
Building mpi_gnu_linux_64
mpifort -m64 -O2 -ffpe-summary=none -cpp -DGITHASH_PP=\"FDS6.6.0-1082-gc728915-master\" -DGITDATE_PP=\""Thu Mar 15 16:33:47 2018 +0100\"" -DBUILDDATE_PP=\""Mar 15, 2018  16:56:09\"" -DCOMPVER_PP=\"unknown\" -DWITH_MKL -I/opt/intel/compilers_and_libraries_2017.5.239/linux/mkl/include -fopenmp -o fds_mpi_gnu_linux_64 prec.o cons.o devc.o data.o type.o mesh.o func.o smvv.o irad.o turb.o soot.o ieva.o pois.o scrc.o radi.o evac.o gsmv.o geom.o part.o vege.o ctrl.o samr.o dump.o hvac.o mass.o read.o wall.o fire.o divg.o velo.o pres.o init.o main.o -Wl,--start-group /opt/intel/compilers_and_libraries_2017.5.239/linux/mkl/lib/intel64/libmkl_gf_lp64.a /opt/intel/compilers_and_libraries_2017.5.239/linux/mkl/lib/intel64/libmkl_gnu_thread.a /opt/intel/compilers_and_libraries_2017.5.239/linux/mkl/lib/intel64/libmkl_core.a /opt/intel/compilers_and_libraries_2017.5.239/linux/mkl/lib/intel64/libmkl_blacs_openmpi_lp64.a -Wl,--end-group -lgomp -lpthread -lm -ldl

``` 




Run FDS case with GLMAT solver

```bash
$ /vagrant/Build/mpi_gnu_linux_64/fds_mpi_gnu_linux_64 /vagrant/Verification/Pressure_Solver/simple_glmat.fds 
  
 Using GLMAT as pressure solver. List of H unknown numbers per proc:
 MYID=       0, NUNKH_LOCAL=    2100

 Fire Dynamics Simulator

 Current Date     : March 15, 2018  16:57:47
 Revision         : FDS6.6.0-1079-g91c42d3-master
 Revision Date    : Wed Feb 14 16:02:36 2018 +0100
 Compiler         : unknown
 Compilation Date : Mar 15, 2018  14:43:09

 MPI Enabled;    Number of MPI Processes:       1
 OpenMP Enabled; Number of OpenMP Threads:      4

 MPI version: 3.1
 MPI library version: Open MPI v3.0.0, package: Open MPI root@vagrant Distribution, ident: 3.0.0, repo rev: v3.0.0, Sep 12, 2017

 Job TITLE        : Compart_Test
 Job ID string    : glmat

 Time Step:      1, Simulation Time:      0.16 s
 Time Step:      2, Simulation Time:      0.32 s
 Time Step:      3, Simulation Time:      0.47 s
 Time Step:      4, Simulation Time:      0.63 s
 Time Step:      5, Simulation Time:      0.79 s
```



 









