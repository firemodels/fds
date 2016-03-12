#!/bin/bash
WHOAMI=`whoami`
if [ ! "$WHOAMI" == "firebot" ]; then
  echo this script must be run as the user firebot
  echo script aborted
  exit
fi

nohup ./run_firebot.sh -u -c -U &
tail -f nohup.out
