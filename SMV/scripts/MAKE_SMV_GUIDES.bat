@echo off

Rem Windows batch file for creating Smokeview Technical Reference guide figures

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

echo Building Smokeview User guide
%svn_drive%
cd %svn_root%\Manuals\SMV_User_Guide
call make_guide

echo Building Smokeview Techncial reference guide
cd %svn_root%\Manuals\SMV_Technical_Reference_Guide
call make_guide

echo Building Smokeview Verification guide
cd %svn_root%\Manuals\SMV_Verification_Guide
call make_guide