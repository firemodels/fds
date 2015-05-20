#!/bin/bash

# expand keyword found in file using newvalue
#   (note that the ~ in %~2 removes surrounding quotes (allowing imbedded blanks in newvalue)

if [ $# -lt 1 ] ; then
  echo "usage: expand_keyword keyword newvalue file"
  echo "       in file, convert all occurrences of \$keyword: .... $ to"
  echo "       \$keyword: newvalue \$"
  exit
fi

keyword=$1
newvalue=$2
file=$3
if ! [ -e $file ]; then
  exit 
fi
sed -e "s/\\\$$keyword:.*\\\$/\\\$$keyword: $newvalue \\\$/g" $file > temp_$$.txt
cp temp_$$.txt $file
rm temp_$$.txt
