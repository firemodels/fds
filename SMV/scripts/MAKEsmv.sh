#!/bin/bash
dir=$1
arg=$2

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd $DIR/../..
GITROOT=`pwd`

cd $GITROOT/SMV/Build/smokeview/$dir
./make_smv.sh $arg
