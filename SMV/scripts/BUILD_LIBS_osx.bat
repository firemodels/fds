@echo off
Title Building 64 bit OSX libraries for smokeview

:: setup environment variables (defining where repository resides etc) 

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

call %envfile%

%svn_drive%
cd %svn_root%\smv\scripts

set scriptdir=FDS-SMV/Utilities/Scripts
set LIBDIR=FDS-SMV/SMV/Build/LIBS

plink %svn_logon% %scriptdir%/ssh_command2.sh %osx_hostname% %LIBDIR%/lib_osx_intel_64 makelibs.sh

echo.
echo build complete
pause
