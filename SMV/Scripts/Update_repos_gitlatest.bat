@echo off

:: batch file used to update Windows, Linux and OSX GIT repos

set envfile="%userprofile%"\fds_smv_env.bat
IF EXIST %envfile% GOTO endif_envexist
echo ***Fatal error.  The environment setup file %envfile% does not exist. 
echo Create a file named %envfile% and use SMV/Scripts/fds_smv_env_template.bat
echo as an example.
echo.
echo Aborting now...
pause>NUL
goto:eof

:endif_envexist

:: location of batch files used to set up Intel compilation environment

call %envfile%

echo.
echo ------------------------------------------------------------------------
echo Updating the Windows repository, %svn_root%, to the latest revision
%svn_drive%
cd %svn_root%
echo Updating the repo:%svn_root%
git remote update
git checkout development
git merge origin/development
git merge firemodels/development

set scriptdir=%linux_svn_root%/SMV/Utilities/Scripts/
set linux_fdsdir=%linux_svn_root%

echo.
echo ------------------------------------------------------------------------
echo Updating the Linux GIT repository, %linux_svn_root%, on %linux_hostname% to the latest revision
plink %linux_logon% %scriptdir%/UPDATE_latest_fds_onhost.csh  %linux_svn_root% %linux_hostname%

echo.
echo ------------------------------------------------------------------------
echo Updating the OSX GIT repository, %linux_svn_root%, on %osx_hostname% to the latest revision
plink %osx_logon% %scriptdir%/UPDATE_latest_fds_onhost.csh  %linux_svn_root% %osx_hostname%

pause
