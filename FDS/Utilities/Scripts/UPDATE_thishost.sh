#!/bin/bash

directory=$1
host=$2

echo Updating the GIT repository $directory on $host to the latest revision
echo
cd ~/$directory
git checkout development
git remote update
git merge origin/development
git merge firemodels/development 
git describe --dirty 
