@echo off

Rem Windows batch file for creating Smokeview guides

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

echo Building Smokeview guides

call %envfile%

%svn_drive%
cd %svn_root%\Manuals\SMV_User_Guide

call make_guide

cd %svn_root%\Manuals\SMV_Technical_Reference_Guide

call make_guide

cd %svn_root%\Manuals\SMV_Verification_Guide

call make_guide