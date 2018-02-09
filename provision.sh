#!/usr/bin/env sh
sudo apt-get update -qq
sudo apt-get install build-essential -y
sudo add-apt-repository ppa:ubuntu-toolchain-r/test -y
sudo apt-get update -qq
sudo apt install gfortran-7 -y
wget https://www.open-mpi.org/software/ompi/v3.0/downloads/openmpi-3.0.0.tar.gz -q
gunzip -c openmpi-3.0.0.tar.gz -q| tar xf -
cd openmpi-3.0.0/
./configure FC=gfortran-7 CC=gcc-7 --prefix=/vagrant/openmpi --enable-static --disable-shared --quiet
make all -j 4 --quiet
make install --quiet
export PATH=/vagrant/openmpi/bin:$PATH
export LD_LIBRARY_PATH=/vagrant/openmpi/lib:$LD_LIBRARY_PATH
cd /vagrant/Build/mpi_gnu_linux_64
./make_fds.sh
