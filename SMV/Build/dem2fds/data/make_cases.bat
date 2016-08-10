@echo off
:: generate terrain with &GEOM
:: set option=-g

:: generate terrain with &OBST
set option=-o

:: generate terrain with &OBST and examine map images
:: set option=-e -o

set dem2fds=..\intel_win_64\dem2fds_win_64.exe
::set dem2fds=dem2fds

echo.
echo demtest
%dem2fds% %option% -n -d %userprofile%\terrain\demtest demtest1.in 
%dem2fds% %option% -n -d %userprofile%\terrain\demtest demtest2.in 

echo.
echo blodget
%dem2fds% %option% -d %userprofile%\terrain\blodget blodget.in 

echo.
echo NIST
%dem2fds% %option% -n -d %userprofile%\terrain\nist nist.in 

echo.
echo tower
%dem2fds% %option% -d %userprofile%\terrain\tower tower.in 

echo.
echo sugarloaf
%dem2fds% %option% -n -d %userprofile%\terrain\sugarloaf sugarloaf.in 

echo.
echo trails
%dem2fds% %option% -d %userprofile%\terrain\trails trails.in 
