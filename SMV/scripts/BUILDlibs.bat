@echo off
set platform=%1

:: batch file for creating libraries on windows, linux or osx

set envfile="%userprofile%"\fds_smv_env.bat
IF EXIST %envfile% GOTO endif_envexist
echo ***Fatal error.  The environment setup file %envfile% does not exist. 
echo Create a file named %envfile% and use SMV/scripts/fds_smv_env_template.bat
echo as an example.
echo.
echo Aborting now...
pause>NUL
goto:eof

:endif_envexist

echo Building Libraries for %platform%
%svn_drive%
set curdir=%CD%

call %envfile%

:: windows

if "%platform%" == "windows" (
title building libraries for windows
cd %svn_root%\SMV\Build\LIBS\intel_win_64
makelibs
goto eof
)

:: osx

if "%platform%" == "osx" (
title building libraries for osx
plink %osx_logon% %linux_svn_root%/SMV/Build/LIBS/intel_osx_64/makelibs.sh
goto eof
)

:: linux

if "%platform%" == "linux" (
title building libraries for linux
plink %linux_logon% %linux_svn_root%/SMV/Build/LIBS/intel_linux_64/makelibs.sh
goto eof
)


:eof
cd %curdir%
pause
