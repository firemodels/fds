#!/bin/bash
dir=$1

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd $DIR/../..
GITROOT=`pwd`

cd $GITROOT/fds/Build/$dir
rm -f *.o *.mod
