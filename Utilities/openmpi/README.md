## OpenMPI Build scripts

This directory contains scripts used to build an OpenMPI version of fds.  Scripts using Intel and Gnu compilers  are found in the `Utilities/openmpi/intel` and `Utilities/openmpi/gnu` directories of this repo. 

The current OSX version of FDS (the Linux and Windows versions of fds are built with the Intel version of MPI) was built with 
[OpenMPI version 3.1.2](https://www.open-mpi.org/software/ompi/v3.1/)
To build OpenMPI, [download OpenMPI 3.1.2](https://download.open-mpi.org/release/open-mpi/v3.1/openmpi-3.1.2.tar.gz) 
( or a different version from [OpenMPI version 3.1.2](https://www.open-mpi.org/software/ompi/v3.1/) ) and follow the instructions/commands below

* `tar xvf openmpi-3.1.2.tar.gz`
* copy CONF_MAKE.sh from intel or gnu repo directory to openmpi-3.1.2 (or the directory created by the tar command)
* edit CONF_MAKE.sh changing --prefix parameter to desired installation directory.  You shouldn't have to change anything else.
* `./CONF_MAKE.sh` to configure and build openmpi
`sudo make install` to complete the installation


