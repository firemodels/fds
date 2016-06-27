@echo off
set dir=%1
set infile=%2

set fulldir=%BASEDIR%\%dir%


cd %fulldir%
echo %infile%
%SMOKEVIEW% -script %infile%_movies.ssf %infile%
