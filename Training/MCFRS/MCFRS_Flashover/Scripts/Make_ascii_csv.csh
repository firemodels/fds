#!/bin/csh -f
#
# Convert a 1D slice file to a csv file.
#
set fds2ascii=~/bin/fds2ascii_linux

cd ../FDS_Input_Files

foreach time (0060 0120 0180)
set low=`expr $time - 5`
set high=`expr $time + 5`
foreach type (00leak 00open)
foreach TC (LL MM RR)
set ans="n"
if ($TC == "LL") set filenum=2
if ($TC == "LL") set ans="z"
if ($TC == "MM") set filenum=3
if ($TC == "RR") set filenum=1
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
