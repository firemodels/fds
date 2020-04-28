#!/bin/bash

# This script saves the openmpi library path and the names of
# all modules loaded when fds was built. 

echo $OPAL_PREFIX > .fdsinfo
echo $LOADEDMODULES | tr ':' ' ' >> .fdsinfo
