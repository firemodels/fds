#!/bin/bash -f
PIDS=`ps -u $USER | grep Xvfb | grep -v grep |  awk '{print $2}'`
for f in $PIDS
do
kill -9 $f
done

