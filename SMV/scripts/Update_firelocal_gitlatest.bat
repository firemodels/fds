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
echo ------------------------------------------------------------------------
echo Updating the Windows repository FIRE-LOCAL to the latest revision
%svn_drive%
cd %userprofile%\FIRE-LOCAL
git remote update
git merge origin/master

set scriptdir=%linux_svn_root%/Utilities/Scripts/

echo.
echo ------------------------------------------------------------------------
echo Updating the Linux GIT repository FIRE-LOCAL on %linux_hostname% to the latest revision
plink %linux_logon% %scriptdir%/UPDATE_repo.sh  FIRE-LOCAL

echo.
echo ------------------------------------------------------------------------
echo Updating the OSX GIT repository FIRE-LOCAL on %osx_hostname% to the latest revision
plink %osx_logon% %scriptdir%/UPDATE_repo.sh  FIRE-LOCAL

pause
