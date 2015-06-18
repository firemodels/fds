#!/bin/bash
CURDIR=`pwd`
repo=FDS-SMVgitclean

cd ~/$repo
git pull
cd $CURDIR
cp $repo/Utilities/Firebot/*.sh .
