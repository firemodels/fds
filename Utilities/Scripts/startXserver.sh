#!/bin/bash -f
lockfile=/tmp/startXlock
GETNEWPORT () 
{
  while [ -e $lockfile ] ; do
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

XVFB=/usr/X11/bin/Xvfb
if [ "`uname`" == "Darwin" ]; then
   XVFB=/usr/X11/bin/Xvfb
fi
GETNEWPORT 
$XVFB :$display_port -fp /usr/share/X11/fonts/misc -screen 0 1280x1024x24 &
export SMV_ID=$!
export DISPLAY=:$display_port
echo $SMV_ID > /tmp/test.out
rm $lockfile
