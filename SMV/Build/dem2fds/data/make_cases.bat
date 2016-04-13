@echo off
:: generate terrain with &GEOM
set option=-e

:: generate terrain with &OBST
::set option=-o

::set dem2fds=..\intel_win_64\dem2fds_win_64.exe
set dem2fds=dem2fds

set case=blodget
%dem2fds% %option% %case% < %case%_elevs.csv > %case%.fds

set case=nist
%dem2fds% %option% %case% < %case%_elevs.csv > %case%.fds

set case=sugarloaf
%dem2fds% %option% %case% < %case%_elevs.csv > %case%.fds

set case=trails
%dem2fds% %option% %case% < %case%_elevs.csv > %case%.fds