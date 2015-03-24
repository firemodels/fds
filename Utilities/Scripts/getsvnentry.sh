#!/bin/bash
DIR=$1
svnnum=`cat ~/FDS-SMV/Validation/$DIR/FDS_Output_Files/*svn.txt 2> /dev/null | sort | head -1`
if [ "$svnnum" != "" ] ; then
  echo $svnnum $DIR
fi
