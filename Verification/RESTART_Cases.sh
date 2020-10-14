#!/bin/bash
OPTION=$1

if [ "$OPTION" == "1" ]; then
  $QFDS -d Miscellaneous restart_test1a.fds
fi
if [ "$OPTION" == "2" ]; then
  $QFDS -d Miscellaneous restart_test1b.fds
fi
