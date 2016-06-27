#!/bin/bash
dir=$1
command=$2
arg=$3

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd $DIR/../..
GITROOT=`pwd`

cd $GITROOT/$dir
./$command $arg
