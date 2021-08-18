#!/bin/bash
for file in *.f90 ;do
  sed 's/[[:blank:]]*$//' $file > xxx.$$
  differ=`diff $file xxx.$$ | wc -l`
  if [ "$differ" != "0" ]; then
    echo $file: trailing white space removed
    mv xxx.$$ $file
  else
     rm -f xxx.$$
  fi
done
