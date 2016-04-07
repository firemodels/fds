@echo off
Title install bundle on a 64 bit OSX system

Rem  Windows batch file to install bundle on a 64 bit OSX system


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

set platform=osx64
plink %osx_logon% %linux_svn_root%/SMV/scripts/run_command.sh Utilities/uploads FDS_%fds_version%-SMV_%smv_version%_%platform%.sh y

echo.
echo installation complete
pause
