#!/bin/bash
DIR=$1
curdir=`pwd`
svnnumdata=/tmp/svnnumdata.$$
ncfile=/tmp/ncfile.$$
cat ~/FDS-SMV/Validation/$DIR/FDS_Output_Files/*svn.txt 2> /dev/null | sort -u > $svnnumdata
svnnum=`head -1 $svnnumdata`
svnnum2=`tail -1 $svnnumdata`
cat $svnnumdata | wc -l > $ncfile 2> /dev/null
nc=`cat $ncfile`
if [ "$svnnum" != "" ] ; then
  cd ~/FDS-SMV/FDS_Source
  svndate=`svn -r $svnnum info | grep Date | awk '{print $4}'`
  if [[ "$nc" -gt 1 ]] ; then
    svnnum=$svnnum-$svnnum2
  fi
  echo "${DIR//_/\_}  & $svndate & $svnnum \\\\ \hline"
  cd $curdir
fi
rm $svnnumdata
rm $ncfile
