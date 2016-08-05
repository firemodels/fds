@echo off
:: generate terrain with &GEOM
:: set option=-g

:: generate terrain with &OBST
set option=-o

set dem2fds=..\intel_win_64\dem2fds_win_64.exe
::set dem2fds=dem2fds

::set case=blodget
::%dem2fds% %option% %case% < %case%_elevs.csv > %case%.fds

::set case=nist
::%dem2fds% %option% %case% < %case%_elevs.csv > %case%.fds

%dem2fds% %option% -d %userprofile%\terrain\tower tower.in 

%dem2fds% %option% -d %userprofile%\terrain\sugarloaf sugarloaf.in 


%dem2fds% %option% -d %userprofile%\terrain\trails trails.in 
