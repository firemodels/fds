#!/bin/bash

directory=$1
host=$2

echo Show revision for the GIT repository $directory on $host
echo
cd ~/$directory
git describe --dirty 
