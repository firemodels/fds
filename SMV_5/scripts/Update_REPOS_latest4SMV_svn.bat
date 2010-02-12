@echo off

Rem Batch file used to update FDS source revision number

set envfile="%homedrive%\%homepath%"\fds_smv_env.bat
IF EXIST %envfile% GOTO endif_envexist
echo ***Fatal error.  The environment setup file %envfile% does not exist. 
echo Create a file named %envfile% and use SMV_5/scripts/fds_smv_env_template.bat
echo as an example.
echo.
echo Aborting now...
pause>NUL
goto:eof

:endif_envexist

Rem location of batch files used to set up Intel compilation environment

call %envfile%

echo.
echo *** Use Smartsvn to update the Windows repository, %svn_root%, to the latest revision

set scriptdir=%linux_svn_root%/Utilities/Scripts/
set linux_fdsdir=%linux_svn_root%

echo.
echo Updating the Linux repository, %linux_svn_root%, to the latest revision
plink %svn_logon% %scriptdir%/UPDATE_latest_fds_onhost.csh  %linux_svn_root% acrux.cfr.nist.gov

echo.
#echo Updating the tiger OSX repository, %linux_svn_root%, to the latest revision
#plink %svn_logon% %scriptdir%/UPDATE_latest_fds_onhost.csh  %linux_svn_root% tiger.cfr.nist.gov
echo Updating the bluesky OSX repository, %linux_svn_root%, to the latest revision
plink %svn_logon% %scriptdir%/UPDATE_latest_fds_onhost.csh  %linux_svn_root% bluesky.cfr.nist.gov

echo.
#echo Updating the devi1 OSX repository, %linux_svn_root%, to the latest revision
#plink %svn_logon% %scriptdir%/UPDATE_latest_fds_onhost.csh  %linux_svn_root% devi1.nist.gov
pause
