@echo off

Rem Windows batch file to upload Smokeview test files to
Rem the download site.  This script assume that the Windows
Rem batch file, MAKEtest.bat, has already been run.

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
cd %svn_root%\SMV\uploads

Rem ----------------------------------------------------------
Rem should not need to edit any lines below

set version=%smv_version%
set platform=win64
set exe=smv_%version%_%platform%.exe

echo Running Smokeview installer: %exe%
pause
call  %exe%


pause
