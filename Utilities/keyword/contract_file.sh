#!/bin/bash

# contract keywords Revision, RevisionDate and CompileDate in file

if [ $# -lt 1 ]; then
  echo usage: contract_file file
  echo        contract all occurrences of the keywords Revision, 
  echo        RevisionDate and CompileDate in file
  exit
fi

bindir=$1
file=$2

if ! [ -e $file ] ; then
  exit 
fi

"$bindir/contract_keyword.sh" Revision $file
"$bindir/contract_keyword.sh" RevisionDate $file
"$bindir/contract_keyword.sh" CompileDate $file
