@echo off
set dir=%1
set infile=%2

set fulldir=%BASEDIR%/%dir%

set in=%infile%
set out=%infile%.err
set stopfile=%infile%.stop

cd %fulldir%
echo **********************
echo %in% started
echo **********************
%CFAST% %in%  
echo **********************
echo %in% completed
echo **********************
