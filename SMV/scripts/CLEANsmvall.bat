@echo off
Title Cleaning Smokeview for 32 bit Windows 

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
echo.
echo cleaning Test Smokeview's on all platforms
call %svn_root%\SMV\scripts\CLEANsmvwin32.bat
call %svn_root%\SMV\scripts\CLEANsmvwin64.bat
call %svn_root%\SMV\scripts\CLEANsmvosx32.bat
call %svn_root%\SMV\scripts\CLEANsmvosx64.bat
call %svn_root%\SMV\scripts\CLEANsmvlinux32.bat
call %svn_root%\SMV\scripts\CLEANsmvlinux64.bat

pause