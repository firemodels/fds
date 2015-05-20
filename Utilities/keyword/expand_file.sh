#!/bin/bash

# expand keywords Revision, RevisionDate and CompileDate in file

if [ $# -lt 1 ] ; then
  echo usage: expand_file dir file
  echo        expand all occurrences of the keywords Revision, 
  echo        RevisionDate and CompileDate in file using properties
  echo        of the directory dir
  exit
fi

bindir=$1
dir=$2
file=$3

fullfile=$dir/$file

if ! [ -e $fullfile ] ; then
  exit 
fi

source "$bindir/get_repo_properties.sh" $dir

"$bindir/expand_keyword.sh" Revision $revision $fullfile
"$bindir/expand_keyword.sh" RevisionDate "$revision_date" $fullfile
"$bindir/expand_keyword.sh" CompileDate "$build_date" $fullfile
