@echo off

Rem Batch file used to update FDS source revision number

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

Rem location of batch files used to set up Intel compilation environment

call %envfile%

echo.
echo Updating the Windows repository, %svn_root%, to the latest revision
%svn_drive%
cd %svn_root%
echo Updating the repo:%svn_root%
git remote update
git checkout development
git merge origin/development
git merge firemodels/development

pause
