#!/bin/bash
pushd `dirname $0` > /dev/null
SCRIPTPATH=`pwd`
popd > /dev/null
echo $SCRIPTPATH
