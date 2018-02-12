#!/usr/bin/env sh
sudo apt-get update -qq
sudo apt-get install build-essential -y
sudo add-apt-repository ppa:ubuntu-toolchain-r/test -y
sudo apt-get update -qq
sudo apt-get install gfortran-7 -y
wget https://www.open-mpi.org/software/ompi/v3.0/downloads/openmpi-3.0.0.tar.gz -q -O /tmp/openmpi-3.0.0.tar.gz
cd /tmp
gunzip -c /tmp/openmpi-3.0.0.tar.gz -q| tar xf -
cd /tmp/openmpi-3.0.0/
./configure FC=gfortran-7 CC=gcc-7 --quiet && make -j 4 --quiet && sudo make install --quiet
#ldconfig -n /vagrant/openmpi/lib
#export PATH=/vagrant/openmpi/bin:$PATH
#cd /vagrant/Build/mpi_gnu_linux_64
#./make_fds.sh
