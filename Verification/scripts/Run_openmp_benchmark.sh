#!/bin/bash
JOBPREFIX=BB_
OUTPUT=/tmp/openmp.out.$$

#---------------------------------------------
#                   wait_cases_release_end
#---------------------------------------------

wait_cases_release_end()
{
   # Scans qstat and waits for cases to end
     while [[ `qstat -a | awk '{print $2 $4}' | grep $(whoami) | grep $JOBPREFIX` != '' ]]; do
        JOBS_REMAINING=`qstat -a | awk '{print $2 $4}' | grep $(whoami) | grep $JOBPREFIX | wc -l`
        echo "Waiting for ${JOBS_REMAINING} openmp benchmark cases to complete."
        sleep 30
     done
}


# run cases

./Run_FDS_Cases.sh -b -j $JOBPREFIX
wait_case_release_end

# copy results to a file
echo `hostname`
echo 64x64x64>$OUTPUT
var=0
for case in a b c d e f g h
do
let "var=var+1"
CSV=../Timing_Benchmarks/openmp_test64${case}_devc.csv
time=`tail -1 $CSV | awk '{print $2}'`
echo $var $time>>$OUTPUT
done
echo "">>$OUTPUT
echo 128x128x128>>$OUTPUT
var=0
for case in a b c d e f g h
do
let "var=var+1"
CSV=../Timing_Benchmarks/openmp_test128${case}_devc.csv
time=`tail -1 $CSV | awk '{print $2}'`
echo $var $time>>$OUTPUT

# output or email results
done
if [ "x$EMAIL" == "x" ]; then
  cat $OUTPUT
else
  cat $OUTPUT | mail -s "openmp benchmark results" $MAILTO
fi
rm $OUTPUT
