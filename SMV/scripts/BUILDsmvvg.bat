@echo off

REM Windows batch file to upload Smokeview test files to
REM the download site.  This script assume that the Windows
REM batch file, MAKEtest.bat, has already been run.

set platform=linux32

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

%svn_drive%

set CURDIR=%CD%

cd %svn_root%\Manuals\SMV_Verification_Guide\
call make_guide

cd %CURDIR%
