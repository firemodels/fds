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
  echo "setting up graphics environment"
  XVFB=Xvfb
  GETNEWPORT 
  $XVFB :$display_port -fp /usr/share/X11/fonts/misc -screen 0 1280x1024x24 &
  export SMV_ID=$!
  export DISPLAY=:$display_port

# Wait for Xvfb
  MAX_ATTEMPTS=240 # About 2 minutes
  COUNT=0
  echo -n "Waiting for Xvfb to be ready..."
  while ! xdpyinfo -display "${DISPLAY}" >/dev/null 2>&1; do
    echo -n "."
    sleep 0.50s
    COUNT=$(( COUNT + 1 ))
    if [ "${COUNT}" -ge "${MAX_ATTEMPTS}" ]; then
      echo "***error: Xserver not initialized on ${DISPLAY}"
      exit 1
    fi
  done
  echo "  Done - Xvfb is ready on $DISPLAY"
  rm -f $lockfile
fi


