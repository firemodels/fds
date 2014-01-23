@echo off

Rem Windows batch file for building Smokeview Verification guide

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

echo Building Smokeview Verification guide

call %envfile%

%svn_drive%
cd %svn_root%\Verification\scripts\

Make_SMV_pictures