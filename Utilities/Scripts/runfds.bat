@echo off

:: $Date$ 
:: $Revision$
:: $Author$

:: set number of openmp threads

set OMP_NUM_THREADS=1

call %SVNROOT%\fds\Utilities\Scripts\getopts.bat %*

set fulldir=%BASEDIR%/%dir%

set in=%infile%.fds
set out=%infile%.err
set stopfile=%infile%.stop

:: test existence of %FDS%

:: test existence of %fulldir%

:: test existence of FDS input file %fulldir%/%in%

cd %fulldir%
echo %in% started

if exist %stopfile% (
   erase %stopfile%
)
if "%rundebug%" == "1" (
   echo 2 > %stopfile%
)

%FDS% %in%  2> %out%

cd %BASEDIR%
