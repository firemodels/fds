@echo off

call %SVNROOT%\Utilities\Scripts\getopts.bat %*

set fulldir=%BASEDIR%/%dir%

cd %fulldir%
echo %infile%
%SMOKEVIEW% -runscript %infile%
