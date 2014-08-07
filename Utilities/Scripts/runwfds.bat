@echo off

call %SVNROOT%\Utilities\Scripts\getopts.bat %*

set fulldir=%BASEDIR%/%dir%

set in=%infile%.fds
set out=%infile%.err
set stopfile=%infile%.stop

Rem test existence of %FDS%

Rem test existence of %fulldir%

Rem test existence of FDS input file %fulldir%/%in%

Rem if STOPFDS=1 then create %fulldir%/%stopfile% and exit

Rem erase %fulldir%\%stopfile%

cd %fulldir%

IF NOT EXIST %WFDSEXE% GOTO notexist
echo %in% started

%WFDS% %in%  > %out%
:notexist

cd %BASEDIR%
