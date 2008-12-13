#!/bin/csh -f
set f2a=~/bin/fds2ascii_linux

foreach n (01 02 03 04 05 06 07)
$f2a < MFRI_Training_Tower_$n\_min_0060.fds2ascii
$f2a < MFRI_Training_Tower_$n\_min_0120.fds2ascii
$f2a < MFRI_Training_Tower_$n\_min_0180.fds2ascii
$f2a < MFRI_Training_Tower_$n\_avg_0060.fds2ascii
$f2a < MFRI_Training_Tower_$n\_avg_0120.fds2ascii
$f2a < MFRI_Training_Tower_$n\_avg_0180.fds2ascii
$f2a < MFRI_Training_Tower_$n\_max_0060.fds2ascii
$f2a < MFRI_Training_Tower_$n\_max_0120.fds2ascii
$f2a < MFRI_Training_Tower_$n\_max_0180.fds2ascii
end
$f2a < MFRI_Training_Tower_01_Savg.fds2ascii
$f2a < MFRI_Training_Tower_01_Lavg.fds2ascii
