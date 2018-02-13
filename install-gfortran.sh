#!/usr/bin/env sh
sudo apt-get update -qq
sudo apt-get install build-essential -y
sudo add-apt-repository ppa:ubuntu-toolchain-r/test -y
sudo apt-get update -qq
sudo apt install gfortran-7 -y