#!/bin/bash

directory=$1
host=$2

cd ~/$directory
git checkout master
git remote update
git merge origin/master
git merge firemodels/master
git describe --dirty --long
