@echo off
:: generate terrain with &GEOM
:: set option=-geom

set dem2fds=..\intel_win_64\dem2fds_win_64.exe
::set dem2fds=dem2fds

%dem2fds% %option% -show -nobuffer -elevdir %userprofile%\terrain\N40W078 -dir %userprofile%\terrain\nist nist.in 
