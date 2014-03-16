#!/bin/bash

echo shutting down graphics environment
sleep 8
if [ "`uname`" == "Darwin" ]; then
  PIDS=`ps -u $USER | grep Xvfb | grep -v grep |  awk '{print $2}'`
  for f in $PIDS
  do
    kill -9 $f
  done
else
  kill $SMV_ID
fi

