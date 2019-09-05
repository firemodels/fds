## OpenMPI Build scripts

This directory contains scripts used to build OpenMPI for the OSX version of fds. 

The current version of FDS was built with 
[OpenMPI version 3.1.2](https://www.open-mpi.org/software/ompi/v3.1/)

To build OpenMPI 3.1.2 

* To build OpenMPI 3.1.2 [download OpenMPI 3.1.2](https://download.open-mpi.org/release/open-mpi/v3.1/openmpi-3.1.2.tar.gz) 
or a different version from [OpenMPI version 3.1.2](https://www.open-mpi.org/software/ompi/v3.1/)
* type: 

`mkdir BUILD`

`cd BUILD`

`tar xvf openmpi-3.1.2.tar.gz`

* copy CONF_MAKE.sh from intel or gnu build directory to BUILD directory

* edit CONF_MAKE.sh changing --prefix parameter to desired installation directory

* type: `./CONF_MAKE.sh` to configure and build openmpi

* type: `sudo make install` to complete the installation


