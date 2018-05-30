@echo off
set CURDIR=%CD%
cd ..\..\..\Build\impi_intel_win_64
git clean -dxf
echo
echo ********** build FDS
echo
call make_fds bot
cd %CURDIR%
