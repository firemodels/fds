#!/bin/bash

CURDIR=`pwd`
SCRIPTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
GITROOT=$SCRIPTDIR/../..
cd $GITROOT
GITROOT=`pwd`

cd $GITROOT/FDS_Source
git clean -dxf
cd $GITROOT/FDS_Compilation
git clean -dxf

cd $GITROOT/SMV/source
git clean -dxf
cd $GITROOT/SMV/Build
git clean -dxf

cd $CURDIR
