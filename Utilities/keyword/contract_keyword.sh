#!/bin/bash

# contract keyword located in file

if [ $# -lt 1 ] ; then
  echo usage: contract_keyword keyword file
  echo "       change all occurences of \$keyword: .... \$ to"
  echo "       \$keyword: unknown \$ in file"
  exit
fi

keyword=$1
file=$2

if ! [ -e $file ]; then
  exit 
fi

sed -e "s/\\\$$keyword:.*\\\$/\\\$$keyword: unknown \\\$/g" $file > temp_$$.txt
cp temp_$$.txt $file 
rm temp_$$.txt
