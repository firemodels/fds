#!/bin/bash
# *** this script will not run if there are untracked or modified 
#     files in your repository

# This script uses commands found at:
# https://help.github.com/articles/dealing-with-line-endings

CURDIR=`pwd`
cd ../..

# check for untracked files

nuntracked=`git clean -dxfn | wc -l`
if [ $nuntracked -gt 0 ]; then
  echo ""
  echo "***Warning:  This repo has $nuntracked untracked files."
  echo Clean the repo before proceeding.
  echo 1.  cd to the repo root
  echo 2.  type: git clean -dxf
  exit
fi

# check for modified files

IS_DIRTY=`git describe --long --dirty | grep dirty | wc -l`
if [ "$IS_DIRTY" == "1" ]; then
  echo ""
  echo "***Warning: This repo has uncommitted changes."
  echo "Commit or revert these changes before proceeding."
  echo "Type: git status -uno"
  echo "to see which files have been changed"
  exit
fi

# failsafe, should't get here if repo has untracked files
# make sure repo is clean (otherwise untracked files will get committed)
while true; do
  read -p "Proceed with line ending update? (yes/no) " yn
  case $yn in
    [Yy]* ) break;;
    [Nn]* ) echo "Update cancelled";exit;;
    * ) echo "Please answer yes or no.";;
  esac
done

git clean -dxf

git add . -u
git commit -m "saving files before refreshing line endings"
git rm --cached -r .
git reset --hard
git add .
git commit -m "normalize all the line endings"
cd $CURDIR
