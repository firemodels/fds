#!/bin/bash
DIR=$1
gitrevisions=/tmp/gitrevisions.$$
ncfile=/tmp/ncfile.$$
cat $FDSSMV/Validation/$DIR/FDS_Output_Files/*git.txt 2> /dev/null | sort -u > $gitrevisions
gitrev=`head -1 $gitrevisions`
gitrev2=`tail -1 $gitrevisions`
cat $gitrevisions | wc -l > $ncfile 2> /dev/null
nc=`cat $ncfile`
if [ "$gitrev" != "" ] ; then
  gitdate=`echo $gitrev | awk -F - '{print $3}' | sed 's/^.\{1\}//' | git show -s --format=%cD  | head -1 | awk '{print $3,$2",",$4}'`
  if [[ "$nc" -gt 1 ]] ; then
    gitrev=$gitrev-$gitrev2
  fi
  echo "${DIR//_/\_}  & $gitdate & $gitrev \\\\ \hline"
fi
rm $gitrevisions
rm $ncfile
