#!/bin/bash
dir=`pwd`
target=${dir##*/}

if [ "$OPAL_PREFIX" != "" ]; then echo $OPAL_PREFIX > .fdsmpi ; fi
echo Building $target
make -j4 VPATH="../../Source" -f ../makefile $target
