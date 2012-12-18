#!/bin/bash -f

OS=`uname`
if [ "$OS" != "Darwin" ]; then
  disp=`id -u`
  Xvfb :$disp -fp /usr/share/X11/fonts/misc -screen 0 1280x1024x24 &
  export SMV_ID=$!
  export DISPLAY=:$disp
  sleep 8
fi
