#!/bin/bash -f
GETNEWPORT () 
{
  display_port=`id -u`
  nmatches=`ps a -e | grep Xvfb | grep $display_port | grep -v grep | wc | awk '{print $1}'`
  echo display_port=$display_port nmatches=$nmatches 
  while [ $nmatches -ne 0 ] ; do
    sleep 1
    display_port=`expr $display_port + 1`
    nmatches=`ps a -e | grep Xvfb | grep $display_port | grep -v grep | wc | awk '{print $1}'`
    echo display_port=$display_port nmatches=$nmatches 
  done
}

OS=`uname`
if [ "$OS" != "Darwin" ]; then
  GETNEWPORT 
  Xvfb :$display_port -fp /usr/share/X11/fonts/misc -screen 0 1280x1024x24 &
  export SMV_ID=$!
  export DISPLAY=:$display_port
  sleep 8
fi
