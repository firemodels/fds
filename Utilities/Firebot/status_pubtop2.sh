#!/bin/bash
CURDIR=`pwd`
cpufrom=~/.firebot/fds_times.csv
sort -n -k 1 -t , $cpufrom | tail -30 | awk -F ',' '{ printf("[%s,%s],\n",$1,$2) }'
cd $CURDIR
