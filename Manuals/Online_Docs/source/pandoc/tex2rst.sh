#!/bin/bash

echo "starting tex2rst ..."

export PANDOC_DIR=`pwd`

cd ../../../FDS_User_Guide/

echo "building reStructuredText ..."
pandoc FDS_User_Guide.tex -o $PANDOC_DIR/test.rst

echo "fixing relative paths for images ..."
sed -i "s/ SCRIPT_FIGURES/ ..\/..\/..\/FDS_User_Guide\/SCRIPT_FIGURES/" $PANDOC_DIR/test.rst
sed -i "s/ FIGURES/ ..\/..\/..\/FDS_User_Guide\/FIGURES/" $PANDOC_DIR/test.rst

echo "finished!"
