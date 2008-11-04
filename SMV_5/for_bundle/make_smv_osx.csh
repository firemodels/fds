#!/bin/csh -f
cd FDS-SMV/MACtiger2/sv5p0
make clean >& /dev/null
make >& ../../bin/make_osx.out
