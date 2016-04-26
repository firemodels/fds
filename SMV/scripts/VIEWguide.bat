@echo off
set guide=%1

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

echo Creating figures for the Smokeview User's guide

call %envfile%

%svn_drive%
cd %svn_root%\Manuals\%guide%

start acrobat %guide%.pdf
