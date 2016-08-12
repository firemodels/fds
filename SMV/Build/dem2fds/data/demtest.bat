@echo off
:: generate terrain with &GEOM
:: set option=-geom

:: generate terrain with &OBST
set option=-obst

set dem2fds=..\intel_win_64\dem2fds_win_64.exe
::set dem2fds=dem2fds

%dem2fds% %option% -nobuffer -show -dir %userprofile%\terrain\demtest demtest1.in 
%dem2fds% %option% -nobuffer -show -dir %userprofile%\terrain\demtest demtest2.in 

