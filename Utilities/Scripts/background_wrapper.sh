#!/bin/bash
script=$1
if [ -e $script ] ; then
  $script
  rm $script
else
  echo $script does not exist
fi
