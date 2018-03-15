#!/usr/bin/env sh

# install build-essential
sudo apt-get install build-essential -qy

# install gfortran 7.X
sudo add-apt-repository ppa:ubuntu-toolchain-r/test -y
sudo apt-get update -qq
sudo apt-get install gfortran-7 git -qy

# compile Open MPI 3.0.0
wget -nv -N https://www.open-mpi.org/software/ompi/v3.0/downloads/openmpi-3.0.0.tar.gz -q
gunzip -c openmpi-3.0.0.tar.gz -q | tar -xf -
cd openmpi-3.0.0
./configure FC=gfortran-7 CC=gcc-7 --prefix=/opt/openmpi-3.0.0 --enable-mpirun-prefix-by-default --enable-static --disable-shared --quiet && make -j 4 --quiet && sudo make install --quiet

# install Intel MKL from apt repository
wget -nv -N https://apt.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS-2019.PUB
sudo apt-key add GPG-PUB-KEY-INTEL-SW-PRODUCTS-2019.PUB -q
sudo sh -c 'echo deb https://apt.repos.intel.com/mkl all main > /etc/apt/sources.list.d/intel-mkl.list'
sudo apt-get update -qq
sudo apt-get install intel-mkl-64bit-2017.4-061 -qy

# compile FDS with GNU/OMPI toolchain with MKL environment for Intel 64 architecture
source /opt/intel/mkl/bin/mklvars.sh intel64
export PATH=/opt/openmpi-3.0.0/bin:$PATH
export LD_LIBRARY_PATH=/opt/openmpi-3.0.0/lib:$LD_LIBRARY_PATH
cd /vagrant/Build/mpi_gnu_linux_64
./make_fds.sh
