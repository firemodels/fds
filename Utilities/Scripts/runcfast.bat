@echo off

call %SVNROOT%\FDS\Utilities\Scripts\getopts.bat %*

set fulldir=%BASEDIR%/%dir%

set in=%infile%
set out=%infile%.err
set stopfile=%infile%.stop

cd %fulldir%
echo %in% started
%CFAST% %in%  
