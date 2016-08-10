@echo off
:: generate terrain with &GEOM
:: set option=-g

:: generate terrain with &OBST
set option=-o

set dem2fds=..\intel_win_64\dem2fds_win_64.exe
::set dem2fds=dem2fds

echo demtest
%dem2fds% %option% -n -e -d %userprofile%\terrain\demtest demtest1.in 
%dem2fds% %option% -n -e -d %userprofile%\terrain\demtest demtest2.in 

