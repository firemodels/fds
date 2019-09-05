## OpenMPI Build scripts

This directory contains scripts used to build OpenMPI for the OSX version of fds.  Scripts using Intel and Gnu compilers  are found in the `Utilities/openmpi/intel` and `Utilities/openmpi/gnu` directories of this repo. 

The current version of FDS was built with 
[OpenMPI version 3.1.2](https://www.open-mpi.org/software/ompi/v3.1/)
To build OpenMPI, [download OpenMPI 3.1.2](https://download.open-mpi.org/release/open-mpi/v3.1/openmpi-3.1.2.tar.gz) 
or a different version from [OpenMPI version 3.1.2](https://www.open-mpi.org/software/ompi/v3.1/) and follow the instructions/commands below

* `mkdir BUILD`
* `cd BUILD`
* copy openmpi-3.1.2.tar.gz into BUILD
* `tar xvf openmpi-3.1.2.tar.gz`
* copy CONF_MAKE.sh from intel or gnu repo directory to BUILD
* edit CONF_MAKE.sh changing --prefix parameter to desired installation directory.  You shouldn't have to change anything else.
* `./CONF_MAKE.sh` to configure and build openmpi
`sudo make install` to complete the installation


