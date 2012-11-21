REM Script to get figures from remote repository
REM uses pscp from the putty distribution

REM change USERHOST and SVNROOT to match your system

set USERHOST=gforney@acrux.cfr.nist.gov
set SVNROOT=FDS-SMV

@echo off
echo getting figures from Linux cluster
erase scriptfigures\*.png
pscp %USERHOST%:%SVNROOT%/Manuals/SMV_5_User_Guide/scriptfigures/*.png scriptfigures\.