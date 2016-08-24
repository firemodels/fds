#!/bin/bash

directory=$1
branch=$2
host=$3

echo
echo Updating the GIT repository: $directory, branch: $branch on host: $host to the latest revision
echo
cd ~/$directory
git checkout $branch
git remote update
git merge origin/$branch
git describe --dirty 
