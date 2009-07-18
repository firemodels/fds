#!/bin/csh -f
#
# Convert a 1D slice file to a csv file.
#
set fds2ascii=~/bin/fds2ascii_linux

cd ../FDS_Input_Files

foreach time (0060 0120 0180)
foreach type (00leak 00open 00case3 00case1)
#foreach time (0060 0120 0180 0240 0300 0360 0420 0480 0540 0600 0660 0720 0780 0840 0900 0960 1020 1080 1140 1200 1260 1320 1380 1440 1500 1560)
#foreach type (00case3m2)
set low=`expr $time - 5`
set high=`expr $time + 5`
foreach TC (LL MM RR)
set ans="n"
if ($TC == "LL") set filenum=1
if ($TC == "LL") set ans="z"
if ($TC == "MM") set filenum=47
if ($TC == "RR") set filenum=48
set fdsin=MCFRS_Flashover\_$type
set csvfile=$fdsin\_$time\_$TC\_f2a.csv
cat << EOF  | $fds2ascii
$fdsin
2
1
$ans
$low $high 0.8
1
$filenum
$csvfile
EOF
end
end
end
