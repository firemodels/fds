#!/bin/csh -f
#
# Convert a 1D slice file to a csv file.
#
set fds2ascii=~/bin/fds2ascii_linux

cd ../FDS_Input_Files

foreach time (0060 0120 0180)
set low=`expr $time - 5`
set high=`expr $time + 5`
foreach n (01 02 03 04 05 06 07)
foreach type (min avg max Savg Lavg)
set fdsin=MFRI_Training_Tower_$n\_$type
set csvfile=$fdsin\_$time\_f2a.csv
set f2a=f2a_$time\_$n\_$type
if ($n != "01" && $type == "Savg" ) goto endtype
if ($n != "01" && $type == "Lavg" ) goto endtype
cat << EOF  | $fds2ascii
$fdsin
2
1
n
$low $high
1
31
$csvfile
EOF
endtype:
end
end
end
