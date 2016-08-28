@echo off
set svn_drive=d:

set SVNROOT=%CD%\..\..\
set WUIDIR=%SVNROOT%\FDS\Verification\Wui

set FDS=%SVNROOT%\FDS\Build\intel_win_64\fds_win_64

echo %FDS%
cd %WUIDIR%
%FDS% onetree_surf_1mesh.fds
smokezip -part2iso onetree_surf_1mesh
pause
