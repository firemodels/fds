 @echo off

:: batch file used to update Windows, Linux and OSX GIT repos

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

:: location of batch files used to set up Intel compilation environment

call %envfile%

echo.
echo ------------------local PC --------------------------------------------
echo Updating the Windows Git repository, %svn_webroot%, to the latest revision

%svn_drive%
cd %svn_webroot%
echo Updating the repo:%svn_webroot%
git remote update
git checkout nist-pages
git merge origin/nist-pages

set scriptdir=%linux_svn_root%/Utilities/Scripts/
set linux_fdsdir=%linux_svn_root%

echo.
echo ------------------ linux: %linux_hostname% --------------------------------------------
plink %linux_logon% %scriptdir%/UPDATE_thishost.sh            %linux_svn_webroot% nist-pages %linux_hostname%

echo.
echo ------------------ linux: %osx_hostname% --------------------------------------------
plink %osx_logon%   %scriptdir%/UPDATE_latest_fds_onhost.csh  %linux_svn_webroot% nist-pages %osx_hostname%

pause
