@echo off
set option=%1
set option2=%2

set dem2fds=..\intel_win_64\dem2fds_win_64.exe
::set dem2fds=dem2fds

%dem2fds% %option% %option2% -nobuffer -dir %userprofile%\terrain\demtest demtest1.in 
%dem2fds% %option% %option2% -nobuffer -dir %userprofile%\terrain\demtest demtest2.in 

%dem2fds% %option% %option2% -dir %userprofile%\terrain\blodget blodget.in 

%dem2fds% %option% %option2% -nobuffer -dir %userprofile%\terrain\nist nist.in 

%dem2fds% %option% %option2% -dir %userprofile%\terrain\tower tower.in 

%dem2fds% %option% %option2% -nobuffer -dir %userprofile%\terrain\sugarloaf sugarloaf.in 

%dem2fds% %option% %option2% -nobuffer -dir %userprofile%\terrain\trails trails.in 
%dem2fds% %option% %option2% -nobuffer -dir %userprofile%\terrain\trails trails2.in 
