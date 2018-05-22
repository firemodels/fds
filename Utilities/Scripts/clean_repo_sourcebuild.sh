#!/bin/bash

CURDIR=`pwd`
SCRIPTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
GITROOT=$SCRIPTDIR/../../..
cd $GITROOT
GITROOT=`pwd`

cd $GITROOT/fds/Source
git clean -dxf

cd $GITROOT/fds/Build
git clean -dxf

cd $GITROOT/fds/Build/Bundle/uploads
git clean -dxf

cd $GITROOT/smv/Source
git clean -dxf

cd $GITROOT/smv/Build
git clean -dxf

cd $GITROOT/smv/Build/Bundle/uploads
git clean -dxf

cd $CURDIR
