@echo off
set fdsrun=..\..\FDS_Compilation\openmp_intel_win_64\fds_openmp_win_64 
call makecase test1
call makecase test2
call makecase test3
call makecase test4
call makecase test5
call makecase test6
call makecase test7
call makecase test8

set OMP_NUM_THREADS=1
%fdsrun% test1.fds

set OMP_NUM_THREADS=2
%fdsrun% test2.fds

set OMP_NUM_THREADS=3
%fdsrun% test3.fds

set OMP_NUM_THREADS=4
%fdsrun% test4.fds

set OMP_NUM_THREADS=5
%fdsrun% test5.fds

set OMP_NUM_THREADS=6
%fdsrun% test6.fds

set OMP_NUM_THREADS=7
%fdsrun% test7.fds

set OMP_NUM_THREADS=8
%fdsrun% test8.fds
