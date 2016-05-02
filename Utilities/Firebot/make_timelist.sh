#!/bin/bash
curdir=`pwd`
cpufrom=~firebot/.firebot/fds_times.csv
historydir=~firebot/.firebot/history
tempfile=/tmp/filelist.$$
# The offset below is computed by substituting
# Jan 1, 2016 5 UTC (12 AM EST) into a web form
# found at:
# http://www.unixtimestamp.com/
BASETIMESTAMP=1451624400

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
