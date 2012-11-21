#!/bin/sh

directory=$1
revision=$2
host=$3

echo Updating directory using
echo directory: $directory
echo  revision: $revision
echo      host: $host
echo
ssh $host \( cd \~/$directory \; svn -r $revision update \)
