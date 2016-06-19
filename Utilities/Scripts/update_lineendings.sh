#!/bin/bash
# *** warning: this script cleans your repo
#     DO NOT RUN if you have any uncommited work you wish to save

# This script uses commands found at:
# https://help.github.com/articles/dealing-with-line-endings

CURDIR=`pwd`
cd ../..
git clean -dxf
git add . -u
git commit -m "saving files before refreshing line endings"
git rm --cached -r .
git reset --hard
git add .
git commit -m "normalize all the line endings"
cd $CURDIR
