#!/bin/bash

directory=$1

echo Updating the GIT repository $directory to the latest revision
echo
cd ~/$directory
git remote update
git merge origin/master
