#!/bin/csh -f
# dummy change
set directory=$1
set      file=$2
set  revision=$3
set      host=$4

echo Updating directory using
echo directory: $directory
echo  revision: $revision
echo      host: $host
echo
ssh $host \( cd \~/$directory \; svn -r $revision update $file \)
