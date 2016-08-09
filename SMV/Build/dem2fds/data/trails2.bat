@echo off
:: generate terrain with &GEOM
:: set option=-g

:: generate terrain with &OBST
set option=-o

set dem2fds=..\intel_win_64\dem2fds_win_64.exe
::set dem2fds=dem2fds

echo.
echo trails
%dem2fds% %option% -e -d %userprofile%\terrain\trails trails2.in 
