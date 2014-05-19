@echo off
set OMP_NUM_THREADS=%1
set input=%2
fds.exe %input%
