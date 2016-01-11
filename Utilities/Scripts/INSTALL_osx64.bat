@echo off
Title Install 64 bit OSX bundle

:: Windows batch file to Install 64 bit OSX bundle

set platform=osx64

set envfile=%userprofile%\fds_smv_env.bat
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

set exe=FDS_%fds_version%-SMV_%smv_version%_%platform%.sh

plink %svn_logon% %linux_svn_root%/SMV/scripts/ssh_command2.sh %osx_hostname% %linux_svn_root%/Utilities/uploads %exe% y
pause
