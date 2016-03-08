#!/bin/bash
curdir=`pwd`
cpufrom=~firebot/.firebot/fds_times.csv
historydir=~firebot/.firebot/history
tempfile=/tmp/filelist.$$
BASETIMESTAMP=1451606400

while getopts 's' OPTION
do
case $OPTION  in
  s)
   cpufrom=~smokebot/.smokebot/smv_times.csv
   historydir=~smokebot/.smokebot/history
   ;;
esac
done
shift $(($OPTIND-1))
for file in $historydir/Git*benchmark*csv 
do
time=`tail -1 $file`
hash=`tail -2 $file | head -1`
gitdate=`git show -s --format=%ct $hash`
gitdate=`echo "scale=5; $gitdate - $BASETIMESTAMP" | bc`
gitdate=`echo "scale=5; $gitdate/86400 " | bc`
echo $gitdate,$time>> $tempfile
done
sort -t ',' -n -k 1 $tempfile
rm $tempfile
