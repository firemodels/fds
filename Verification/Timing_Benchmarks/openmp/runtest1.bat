@echo off

set FDS=%USERPROFILE%\FDS-SMV\FDS_Compilation\openmp_intel_win_64\fds_openmp_win_64

set OMP_NUM_THREADS=1
call %FDS% test1_01.fds

set OMP_NUM_THREADS=2
call %FDS% test1_02.fds

set OMP_NUM_THREADS=4
call %FDS% test1_04.fds

set OMP_NUM_THREADS=8
call %FDS% test1_08.fds

set OMP_NUM_THREADS=12
call %FDS% test1_12.fds

set OMP_NUM_THREADS=16
call %FDS% test1_16.fds