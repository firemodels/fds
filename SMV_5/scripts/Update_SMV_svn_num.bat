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

%svn_drive%
set win_dir=%svn_root%\SMV_5\source\smokeview
cd %win_dir%

echo.
echo Updating the directory %win_dir% on the Windows repository to the SVN revision: %smv_revision%
svn -r %smv_revision% update

set scriptdir=%linux_svn_root%/SMV_5/scripts/
set linux_smvdir=%linux_svn_root%/SMV_5/source/smokeview

echo.
echo Updating the directory %linux_smvdir% on the Linux repository to the SVN revision: %smv_revision%
plink %svn_logon% %scriptdir%/UPDATE_smv_onhost.csh  %linux_smvdir% %smv_revision% acrux.cfr.nist.gov

echo.
echo Updating the directory %linux_smvdir% on the the OSX repository to the SVN revision: %smv_revision%
plink %svn_logon% %scriptdir%/UPDATE_smv_onhost.csh  %linux_smvdir% %smv_revision% tiger.cfr.nist.gov

pause