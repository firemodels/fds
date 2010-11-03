@echo off

Rem Batch file used to update FDS source revision number

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

Rem location of batch files used to set up Intel compilation environment

call %envfile%

%svn_drive%
set win_fdsdir=%svn_root%\FDS_Source
cd %win_fdsdir%

echo.
echo *** Use Smartsvn to update %win_fdsdir% in the Windows repository to SVN revision: %fds_revision%

set scriptdir=%linux_svn_root%/FDS_Compilation/Scripts/
set linux_fdsdir=%linux_svn_root%/FDS_Source

echo.
echo Updating %linux_fdsdir% in the Linux repository on acrux to SVN revision: %fds_revision%
plink %svn_logon% %scriptdir%/UPDATE_fds_onhost.csh  %linux_fdsdir% %fds_revision% acrux.cfr.nist.gov

echo.
echo Updating %linux_fdsdir% in the OSX repository on %OSXHOST% to SVN revision: %fds_revision%
plink %svn_logon% %scriptdir%/UPDATE_fds_onhost.csh  %linux_fdsdir% %fds_revision% %OSXHOST%

pause
