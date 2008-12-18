#!/bin/csh -f
#
#  merge 3 door opening time cases into 1 csv file
#
set pid=$$
set headskip=~/bin/headskip
set catcol=~/bin/catcol
set INDIR=../FDS_Input_Files
set OUTDIR=../FDS_Output_Files

set file1=MFRI_Training_Tower_01_Savg_devc.csv
set file2=MFRI_Training_Tower_01_avg_devc.csv
set file3=MFRI_Training_Tower_01_Lavg_devc.csv
set file1t=aa.$$
set file2t=bb.$$
set file3t=cc.$$
set fileout=$OUTDIR/MFRI_Training_Tower_door_time.csv

cd $INDIR
echo s,C,C,C,C,C,C,C,C,C > $fileout
echo FDS Time,TC03S,TC06S,TC10S,TC03M,TC06M,TC10M,TC03L,TC06L,TC10L>>$fileout
$headskip -s 2 < $file1 | cut -d , -f 1,5,8,12 - >$file1t
$headskip -s 2 < $file2 | cut -d , -f 5,8,12 - >$file2t
$headskip -s 2 < $file2 | cut -d , -f 5,8,12 - >$file3t

$catcol $file1t $file2t $file3t >> $fileout
rm $file1t $file2t $file3t
