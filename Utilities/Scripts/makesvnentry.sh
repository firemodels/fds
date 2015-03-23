#!/bin/bash
DIR=$1
curdir=`pwd`
svnnumdata=/tmp/svnnumdata.$$
cat ~/FDS-SMV/Validation/$DIR/FDS_Output_Files/*svn.txt > $svnnumdata 2> /dev/null
svnnum=`tail -1 $svnnumdata`
if [ "$svnnum" != "" ] ; then
  cd ~/FDS-SMV/FDS_Source
  svndate=`svn -r $svnnum info | grep Date | awk '{print $4}'`
  echo "${DIR//_/\_}  & $svndate & $svnnum \\\\ \hline"
  cd $curdir
fi
rm $svnnumdata
