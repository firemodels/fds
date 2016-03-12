#!/bin/bash
WHOAMI=`whoami`
if [ ! "$WHOAMI" == "firebot" ]; then
  echo this script must be run as the user firebot
  echo script aborted
  exit
fi
echo You are about to run firebot.  
echo "Press any key to proceed or <CTRL> c to abort."
read val

nohup ./run_firebot.sh -u -c -U &
echo firebot started.
echo
echo Type: 
echo tail -f nohup.out 
echo if you wish to see the status of the firebot while it is running.
echo Look in the output directory for files named warnings or errors
echo to see interim firebot results.
