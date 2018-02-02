#!/usr/bin/env sh
sudo apt-get update -qq
sudo apt-get install build-essential -y
sudo add-apt-repository ppa:ubuntu-toolchain-r/test -y
sudo apt-get update -qq
sudo apt install gfortran-7 -y
wget https://www.open-mpi.org/software/ompi/v3.0/downloads/openmpi-3.0.0.tar.gz -q
#consider: tar -xzf
gunzip -c openmpi-3.0.0.tar.gz -q| tar xf -
cd openmpi-3.0.0/
#./configure FC=gfortran-7 CC=gcc-7 --prefix=/opt/openmpi_64 --enable-mpirun-prefix-by-default --enable-mpi-fortran --enable-static --disable-shared --quiet
# consider: build shared libraries
./configure FC=gfortran-7 CC=gcc-7 --prefix=/opt/openmpi_64 --quiet
make -j 4 --quiet
sudo make install --quiet
# why this is necessary is explain in OpenMPI README file
# using absolute path to mpirun is an alternative
export PATH=/opt/openmpi_64/bin:$PATH
export LD_LIBRARY_PATH=/opt/openmpi_64/lib:$LD_LIBRARY_PATH
cd /vagrant/
cd Build/mpi_gnu_linux_64
./make_fds.sh
