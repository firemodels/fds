@echo off

Rem Batch file used to convert FDS wiki to html

set envfile="%userprofile%\fds_smv_env.bat"
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

%svn_drive%
cd %svn_root%\Utilities\Scripts\bundle_setup

set scriptdir=FDS-SMV/Utilities/Scripts
set bundlesetup=%scriptdir%/bundle_setup


echo.
echo Converting the FDS release notes from wiki to html format

plink %svn_logon%  %scriptdir%/ssh_command2.csh %linux_hostname% %scriptdir% CONV_fds_release.sh

pscp %svn_logon%:%bundlesetup%/release_notes.htm FDS_Release_Notes.htm
pscp %svn_logon%:%bundlesetup%/fds_linux.html README_LINUX.html
pscp %svn_logon%:%bundlesetup%/fds_osx.html   README_OSX.html
pause
