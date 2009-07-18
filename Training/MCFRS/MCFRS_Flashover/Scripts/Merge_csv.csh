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
foreach test (00leak 00open 00case1 00case3)
#foreach time (0060 0120 0180 0240 0300 0360 0420 0480 0540 0600 0660 0720 0780 0840 0900 0960 1020 1080 1140 1200 1260 1320 1380 1440 1500 1560)
#foreach test (00case3m2)
echo merging test $test, time $time
set file1=MCFRS_Flashover_$test\_$time\_LL\_f2a.csv
set file2=MCFRS_Flashover_$test\_$time\_MM\_f2a.csv
set file3=MCFRS_Flashover_$test\_$time\_RR\_f2a.csv
set file1t=aa.$$
set file2t=bb.$$
set file3t=cc.$$
set fileout=$OUTDIR/MCFRS_Flashover_$test\_all_$time\_f2a.csv

echo Z,TEMP_LL,TEMP_MM,TEMP_RR>$fileout
$headskip -s 1 <$file1 >$file1t
$headskip -s 1 <$file2 | cut -d , -f 2 - >$file2t
$headskip -s 1 <$file3 | cut -d , -f 2 - >$file3t

$catcol $file1t $file2t $file3t >> $fileout
rm $file1t $file2t $file3t
end
end
