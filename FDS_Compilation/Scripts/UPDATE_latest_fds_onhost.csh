#!/bin/csh -f
# dummy
set directory=$1
set      host=$2

echo SVN updating using
echo directory: $directory
echo      host: $host
echo
ssh $host \( cd \~/$directory \; svn update \)
