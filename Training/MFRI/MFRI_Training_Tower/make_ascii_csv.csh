#!/bin/csh -f
set f2a=~/bin/fds2ascii_linux

$f2a < MFRI_Training_Tower_01_min.fds2ascii
$f2a < MFRI_Training_Tower_02_min.fds2ascii
#$f2a < MFRI_Training_Tower_03_min.fds2ascii
$f2a < MFRI_Training_Tower_04_min.fds2ascii
$f2a < MFRI_Training_Tower_05_min.fds2ascii
$f2a < MFRI_Training_Tower_06_min.fds2ascii
$f2a < MFRI_Training_Tower_07_min.fds2ascii

$f2a < MFRI_Training_Tower_01_avg.fds2ascii
$f2a < MFRI_Training_Tower_02_avg.fds2ascii
$f2a < MFRI_Training_Tower_03_avg.fds2ascii
$f2a < MFRI_Training_Tower_04_avg.fds2ascii
$f2a < MFRI_Training_Tower_05_avg.fds2ascii
$f2a < MFRI_Training_Tower_06_avg.fds2ascii
$f2a < MFRI_Training_Tower_07_avg.fds2ascii

#$f2a < MFRI_Training_Tower_01_max.fds2ascii
$f2a < MFRI_Training_Tower_02_max.fds2ascii
$f2a < MFRI_Training_Tower_03_max.fds2ascii
$f2a < MFRI_Training_Tower_04_max.fds2ascii
$f2a < MFRI_Training_Tower_05_max.fds2ascii
$f2a < MFRI_Training_Tower_06_max.fds2ascii
$f2a < MFRI_Training_Tower_07_max.fds2ascii
