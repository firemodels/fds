#!/bin/bash

echo "starting tex2rst ..."

export PANDOC_DIR=`pwd`

cd ../../../FDS_User_Guide/
pandoc FDS_User_Guide.tex -o $PANDOC_DIR/test.rst

echo "finished!"
