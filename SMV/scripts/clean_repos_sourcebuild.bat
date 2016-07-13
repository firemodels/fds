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
echo ------------------------------------------------------------------------
echo Cleaning source and build directories in the Windows repository %svn_root%
%svn_drive%
cd %svn_root%\SMV\Build
git clean -dxf
cd %svn_root%\SMV\Source
git clean -dxf

set scriptdir=%linux_svn_root%/Utilities/Scripts/

echo.
echo ------------------------------------------------------------------------
echo Cleaning source and build directories in the Linux repository %linux_svn_root%, on %linux_hostname%
plink %linux_logon% %scriptdir%/clean_repo_sourcebuild.sh  %linux_svn_root% %linux_hostname%

echo.
echo ------------------------------------------------------------------------
echo Cleaning source and build directories in the OSX repository %linux_svn_root%, on %osx_hostname%
plink %osx_logon% %scriptdir%/clean_repo_sourcebuild.sh

pause
