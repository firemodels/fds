I. Build OpenMPI library

0. cd into OPENMPI

1. download the openmpi-1.8.4.tar.gz from
http://www.open-mpi.org/software/ompi/v1.8/

2. untar openmpi-1.8.4.tar.gz and copy CONFIGURE.sh, MAKE.sh and MAKEINSTALL.sh into the top level of what you just untar'd

3. examine CONFIGURE.sh, MAKE.sh, MAKEINSTALL.sh and change compiler location if necessary.

4. In CONFIGURE.sh, the only other setting you should need to change is --prefix, the location and name of the installed openmpi
library

5. run CONFIGURE.sh, MAKE.sh as a normal user

6. run MAKEINSTALL.sh as root if the install directory is owned by root

II.  Examine .bashrc and change compiler and openmpi location if necessary

