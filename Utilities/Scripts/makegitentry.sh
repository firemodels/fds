#!/bin/bash
DIR=$1
gitrevisions=/tmp/gitrevisions.$$
cat $FDSSMV/Validation/$DIR/FDS_Output_Files/*git.txt 2> /dev/null | sort -u > $gitrevisions
gitrev=`head -1 $gitrevisions`
if [ "$gitrev" != "" ] ; then
  gitrevshort=`echo $gitrev | awk -F - '{print $3}' | sed 's/^.\{1\}//'`
  gitdate=`git show -s --format=%aD $gitrevshort | head -1 | awk '{print $3,$2",",$4}'`
  gitdate2=`git show -s --format=%at $gitrevshort | head -1 | awk '{print $1}'`
  echo "${DIR//_/\_}  & $gitdate & $gitrev & $gitdate2\\\\ \hline"
fi
rm $gitrevisions
