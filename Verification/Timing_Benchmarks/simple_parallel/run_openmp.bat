@echo off
set fdsrun=..\..\..\FDS_Compilation\openmp_intel_win_64\fds_openmp_win_64 
call makecase bencho1
call makecase bencho2
call makecase bencho3
call makecase bencho4
call makecase bencho5
call makecase bencho6
call makecase bencho7
call makecase bencho8

set OMP_NUM_THREADS=1
%fdsrun% bencho1.fds

set OMP_NUM_THREADS=2
%fdsrun% bencho2.fds

set OMP_NUM_THREADS=3
%fdsrun% bencho3.fds

set OMP_NUM_THREADS=4
%fdsrun% bencho4.fds

set OMP_NUM_THREADS=5
%fdsrun% bencho5.fds

set OMP_NUM_THREADS=6
%fdsrun% bencho6.fds

set OMP_NUM_THREADS=7
%fdsrun% bencho7.fds

set OMP_NUM_THREADS=8
%fdsrun% bencho8.fds
