I. Build OpenMPI library

0. cd into OPENMPI

1. download openmpi-1.8.4.tar.gz from 
http://www.open-mpi.org/software/ompi/v1.8/

2. untar openmpi-1.8.4.tar.gz and copy CONFIGURE.sh, MAKE.sh and MAKEINSTALL.sh into openmpi-1.8.4 

3. cd into openmpi-1.8.4 and examine CONFIGURE.sh, MAKE.sh, MAKEINSTALL.sh and change compiler location if necessary.

4. In CONFIGURE.sh, the only other setting you should need to change is --prefix, the location and name of the installed openmpi
library

5. run CONFIGURE.sh, MAKE.sh as a normal user

6. run MAKEINSTALL.sh as root if the install directory is owned by root


II.  Build FDS (Fire Dynamics Simulator)

1. Examine .bashrc and change compiler and openmpi location if necessary

2. cd to ~/FDS-SMV/FDS_Compilation/mpi_intel_linux_64ib

3. execute the script make_fds.sh


III. cd back to ~/FDS-SMV/Training/Miscellaneous

1. run the command: ~/FDS-SMV/Utilities/Scripts/qfds.sh -p 120 -n 6 test_mpi_120.fds

2. check the q to see if the job is running

