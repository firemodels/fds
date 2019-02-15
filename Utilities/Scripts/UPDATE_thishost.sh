#!/bin/bash

directory=$1

cd ~/$directory
git checkout master
git remote update
git merge origin/master
git merge firemodels/master
git describe --dirty --long
