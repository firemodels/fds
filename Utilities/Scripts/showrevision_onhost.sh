#!/bin/bash

directory=$1
host=$2

echo Show revision for the GIT repository $directory on $host
echo
ssh -q $host \( cd \~/$directory \; git describe --dirty  \)
