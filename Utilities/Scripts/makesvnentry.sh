#!/bin/bash
DIR=$1
curdir=`pwd`
svnnumdata=/tmp/svnnumdata.$$
svndatedata=/tmp/svndatedata.$$
ncfile=/tmp/ncfile.$$
cat ~/FDS-SMV/Validation/$DIR/FDS_Output_Files/*svn.txt 2> /dev/null | awk '{print $1}' | sort -u > $svnnumdata
cat ~/FDS-SMV/Validation/$DIR/FDS_Output_Files/*svn.txt 2> /dev/null | awk '{print $3}' | sort -u > $svndatedata
svnnum=`head -1 $svnnumdata`
svnnum2=`tail -1 $svnnumdata`
cat $svnnumdata | wc -l > $ncfile 2> /dev/null
nc=`cat $ncfile`
if [ "$svnnum" != "" ] ; then
  svndate=`head -1 $svndatedata`
  if [[ "$nc" -gt 1 ]] ; then
    svnnum=$svnnum-$svnnum2
  fi
  echo "${DIR//_/\_}  & $svndate & $svnnum \\\\ \hline"
  cd $curdir
fi
rm $svnnumdata
rm $svndatedata
rm $ncfile
