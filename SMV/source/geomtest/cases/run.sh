#!/bin/bash
export GEOMTEST=../intel_linux_64/geomtest
OS=`uname`
if [ "$OS" == "Darwin" ]; then
  export GEOMTEST=../intel_osx_64/geomtest
fi

./geom_cases.sh
