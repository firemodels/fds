#!/bin/csh -f
set host=$1
set fromdir=$2
set todir=$3
scp -r $fromdir/* $host\:$todir
