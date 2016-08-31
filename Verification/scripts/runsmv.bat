@echo off

call %SVNROOT%\fds\Utilities\Scripts\getopts.bat %*

set fulldir=%BASEDIR%/%dir%

cd %fulldir%
echo %infile%
%SMOKEVIEW% -runscript %infile%
