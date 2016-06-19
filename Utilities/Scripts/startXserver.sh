#!/bin/bash
lockfile=/tmp/startXlock
GETNEWPORT () 
{
  while [ -e $lockfile ] ; do
    echo waiting for lock file, $lockfile, to clear
    sleep 5
  done
  touch $lockfile
  chmod 777 $lockfile
  display_port=`id -u`
  nmatches=`ps a -e | grep $XVFB | grep $display_port | grep -v grep | wc | awk '{print $1}'`
  while [ $nmatches -ne 0 ] ; do
    display_port=`expr $display_port + 1`
    nmatches=`ps a -e | grep $XVFB | grep $display_port | grep -v grep | wc | awk '{print $1}'`
  done
}

if [ "`uname`" != "Darwin" ]; then
  PAUSE=30
  echo "setting up graphics environment (pausing $PAUSE s)"
  XVFB=Xvfb
  GETNEWPORT 
  $XVFB :$display_port -fp /usr/share/X11/fonts/misc -screen 0 1280x1024x24 &
  export SMV_ID=$!
  export DISPLAY=:$display_port
  rm -f $lockfile
  sleep $PAUSE
fi
