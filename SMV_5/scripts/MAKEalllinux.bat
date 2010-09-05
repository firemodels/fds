@echo off

Rem  Windows batch file to copy smokediffset SVNROOT=~/FDS-SMV

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

call %svn_root%\SMV_5\scripts\MAKEsmdlinux.bat
call %svn_root%\SMV_5\scripts\MAKEsmzlinux.bat
call %svn_root%\SMV_5\scripts\MAKEbglinux.bat
call %svn_root%\SMV_5\scripts\MAKEf2alinux.bat
