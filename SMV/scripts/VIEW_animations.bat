@echo off

Rem Windows batch file for creating Smokeview User guide figures

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

echo Viewing animation web page

call %envfile%

%svn_drive%

start explorer %svn_root%\Manuals\SMV_Animations\index.html
