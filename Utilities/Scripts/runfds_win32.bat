@echo off
set dir=%1
set infile=%2

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
%FDS% %in%  > %out%
