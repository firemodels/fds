@echo off
set option=%1
set option2=%2

set dem2fds=..\intel_win_64\dem2fds_win_64.exe
::set dem2fds=dem2fds

echo.
echo demtest
%dem2fds% %option% %option2% -nobuffer -dir %userprofile%\terrain\demtest demtest1.in 
%dem2fds% %option% %option2% -nobuffer -dir %userprofile%\terrain\demtest demtest2.in 

echo.
echo blodget
%dem2fds% %option% %option2% -dir %userprofile%\terrain\blodget blodget.in 

echo.
echo NIST
%dem2fds% %option% %option2% -nobuffer -dir %userprofile%\terrain\nist nist.in 

echo.
echo tower
%dem2fds% %option% %option2% -dir %userprofile%\terrain\tower tower.in 

echo.
echo sugarloaf
%dem2fds% %option% %option2% -nobuffer -dir %userprofile%\terrain\sugarloaf sugarloaf.in 

echo.
echo trails
%dem2fds% %option% %option2% -dir %userprofile%\terrain\trails trails.in 
