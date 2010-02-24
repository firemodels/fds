#!/bin/csh -f
# dummy change
set directory=$1
set  revision=$2
set      host=$3

echo Updating directory using
echo directory: $directory
echo  revision: $revision
echo      host: $host
echo
ssh $host \( cd \~/$directory \; svn -r $revision update \)
