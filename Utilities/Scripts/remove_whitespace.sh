#!/bin/bash
for file in *.f90 ;do
  sed 's/[[:blank:]]*$//' $file > xxx.$$
  nlines=`diff $file xxx.$$ | wc -l`
  if [ "$nlines" != "0" ]; then
    echo $file: trailing white space removed from $nlines lines
    mv xxx.$$ $file
  fi
done
