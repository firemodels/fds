#!/bin/bash
DIR=$1
curdir=`pwd`
gitrevisions=/tmp/gitrevisions.$$
ncfile=/tmp/ncfile.$$
cat $FDSSMV/Validation/$DIR/FDS_Output_Files/*git.txt 2> /dev/null | sort -u > $gitrevisions
gitrev=`head -1 $gitrevisions`
gitrev2=`tail -1 $gitrevisions`
cat $gitrevisions | wc -l > $ncfile 2> /dev/null
nc=`cat $ncfile`
if [ "$gitrev" != "" ] ; then
  cd $FDSSMV/FDS_Source
  gitdate=`git log .  | head -3 | tail -1 | awk '{print $3,$4",",$6}'`
  if [[ "$nc" -gt 1 ]] ; then
    gitrev=$gitrev-$gitrev2
  fi
  echo "${DIR//_/\_}  & $gitdate & $gitrev \\\\ \hline"
  cd $curdir
fi
rm $gitrevisions
rm $ncfile
