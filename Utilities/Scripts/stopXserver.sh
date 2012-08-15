#!/bin/bash -f

OS=`uname`
if [ "$OS" != "Darwin" ]; then
  sleep 8
  kill $SMV_ID
fi
