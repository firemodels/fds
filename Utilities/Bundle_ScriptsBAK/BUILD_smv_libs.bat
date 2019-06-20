@echo off

:: batch file for creating libraries on windows, linux or osx

set envfile="%userprofile%"\fds_smv_env.bat
IF EXIST %envfile% GOTO endif_envexist
echo ***Fatal error.  The environment setup file %envfile% does not exist. 
echo Create a file named %envfile% and use smv/scripts/fds_smv_env_template.bat
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

title building libraries for windows
cd %svn_root%\smv\Build\LIBS\intel_win_64
make_LIBS

:: osx

title building libraries for osx
plink %osx_logon% %linux_svn_root%/smv/Build/LIBS/intel_osx_64/make_LIBS.sh

:: linux

title building libraries for linux
plink %linux_logon% %linux_svn_root%/smv/Build/LIBS/intel_linux_64/make_LIBS.sh

cd %curdir%
pause
