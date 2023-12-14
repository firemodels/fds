#!/bin/bash
# Command Line Tools for Xcode
short=0
if [[ $(type -t pkgutil) == file ]]; then
if pkgutil --pkgs=com.apple.pkg.CLTools_Executables >/dev/null; then
    long=$(pkgutil --pkg-info=com.apple.pkg.CLTools_Executables | awk '/version:/ {print $2}')
    short="${long:0:2}"
fi
fi
if [[ "$short" == "15" ]]; then
   echo "-ld_classic"
else
   echo " "
fi
