@echo off

Rem Windows batch file to upload Smokeview test files to
Rem the google download site.  This script assume that the Windows
Rem batch file, MAKEtest.bat, has already been run.

set envfile="%homedrive%\%homepath%"\fds_smv_env.bat
IF EXIST %envfile% GOTO endif_envexist
echo ***Fatal error.  The environment setup file %envfile% does not exist. 
echo Create a file named %envfile% and use SMV/scripts/fds_smv_env_template.bat
echo as an example.
echo.
echo Aborting now...
pause>NUL
goto:eof

:endif_envexist

set platform=win32

call %envfile%

echo Running Smokeview installer:  smv_test_%smv_revision%_%platform%.exe
pause

%svn_drive%
cd %svn_root%\smv_5\for_bundle\to_google

call smv_test_%smv_revision%_%platform%.exe
pause
