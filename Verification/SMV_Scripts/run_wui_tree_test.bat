@echo off
set svn_drive=d:

set SVNROOT=%CD%\..\..\
set WUIDIR=%SVNROOT%\Verification\Wui

set FDS=%SVNROOT%\FDS_Compilation\intel_win_32\fds5_win_32

echo %FDS%
cd %WUIDIR%
%FDS% onetree_surf_1mesh.fds
smokezip -part2iso onetree_surf_1mesh
pause
