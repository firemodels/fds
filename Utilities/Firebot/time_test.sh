#!/bin/bash

GET_TIME(){
  echo $(date +"%s")
}
GET_DURATION(){
  time_before=$1
  time_after=$2
  DIFF_TIME=`echo $(($time_after-$time_before))`
  echo "$(($DIFF_TIME / 3600 ))h $((($DIFF_TIME % 3600) / 60))m $(($DIFF_TIME % 60))s"
}

s1=`GET_TIME`
echo time1=$s1
sleep 70
s2=`GET_TIME`
echo time2=$s2
echo duration=`GET_DURATION $s1 $s2`
