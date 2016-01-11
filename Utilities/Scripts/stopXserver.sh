#!/bin/bash
if [ "`uname`" != "Darwin" ]; then
  echo shutting down graphics environment
  sleep 8
  kill $SMV_ID
fi
