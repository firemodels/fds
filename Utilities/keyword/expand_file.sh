#!/bin/bash

# expand keywords Revision, RevisionDate and CompileDate in file

if [ $# -lt 1 ] ; then
  echo usage: expand_file bindir dir file
  echo        expand all occurrences of the keywords Revision, 
  echo        RevisionDate and CompileDate in file using properties
  echo        of the directory dir
  exit
fi

bindir=$1
dir=$2
file=$3

if ! [ -e $file ] ; then
  exit 
fi

source "$bindir/get_repo_properties.sh" $dir

"$bindir/expand_keyword.sh" Revision $revision $file
"$bindir/expand_keyword.sh" RevisionDate "$revision_date" $file
"$bindir/expand_keyword.sh" CompileDate "$build_date" $file
