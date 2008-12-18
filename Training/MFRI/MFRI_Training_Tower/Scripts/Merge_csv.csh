#!/bin/csh -f
#
# merge 3 csv files  for 3 simulation times into one csv file
#
set pid=$$
set headskip=~/bin/headskip
set catcol=~/bin/catcol
set INDIR=../FDS_Input_Files
set OUTDIR=../FDS_Output_Files
cd $INDIR
foreach time (0060 0120 0180)
foreach test (01 02 03 04 05 06 07)
echo merging test $test, time $time
set file1=MFRI_Training_Tower_$test\_min_$time\_f2a.csv
set file2=MFRI_Training_Tower_$test\_avg_$time\_f2a.csv
set file3=MFRI_Training_Tower_$test\_max_$time\_f2a.csv
set file1t=aa.$$
set file2t=bb.$$
set file3t=cc.$$
set fileout=$OUTDIR/MFRI_Training_Tower_$test\_all_$time\_f2a.csv

echo Z,TEMP_MIN,TEMP_AVG,TEMP_MAX>$fileout
$headskip -s 1 <$file1 >$file1t
$headskip -s 1 <$file2 | cut -d , -f 2 - >$file2t
$headskip -s 1 <$file3 | cut -d , -f 2 - >$file3t

$catcol $file1t $file2t $file3t >> $fileout
rm $file1t $file2t $file3t
end
end
