@echo off

Rem  Windows batch file to copy smokezip

Rem setup environment variables (defining where repository resides etc) 

set envfile="%homedrive%\%homepath%"\fds_smv_env.bat
IF EXIST %envfile% GOTO endif_envexist
echo ***Fatal error.  The environment setup file %envfile% does not exist. 
echo Create a file named %envfile% and use SMV_5/scripts/fds_smv_env_template.bat
echo as an example.
echo.
echo Aborting now...
pause>NUL
goto:eof

:endif_envexist

call %envfile%

%svn_drive%
cd %svn_root%\smv_5\scripts

set scriptdir=FDS-SMV/SMV_5/scripts

plink %svn_logon% %scriptdir%/ssh_command.csh fire79 %scriptdir% COPYsmzlinux.csh

echo.
echo copy complete
pause
