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
foreach test (00leak 00open)
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
