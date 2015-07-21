#!/bin/csh -f

set directory=$1
set      host=$2

echo SVN updating using
echo directory: $directory
echo      host: $host
echo
ssh -q $host \( cd \~/$directory \; git checkout development \; git remote update \; git pull \)
