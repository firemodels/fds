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

%svn_drive%
set win_fdsdir=%svn_root%\FDS_Source
cd %win_fdsdir%

echo.
echo *** Use Smartsvn to update %win_fdsdir% in the Windows repository to SVN revision: %fds_revision%

set scriptdir=%linux_svn_root%/FDS_Compilation/Scripts/
set fds_sourcedir=%linux_svn_root%/FDS_Source
set fds_builddir=%linux_svn_root%/FDS_Compilation

echo.
echo Updating %fds_sourcedir% in the Linux repository on %linux_hostname$ to SVN revision: %fds_revision%
plink %svn_logon% %scriptdir%/UPDATE_fds_onhost.csh  %fds_sourcedir% %fds_revision% %linux_hostname%
plink %svn_logon% %scriptdir%/UPDATE_fds_onhost_file.csh  %fds_builddir% makefile %fds_revision% %linux_hostname%

echo.
echo Updating %fds_sourcedir% in the OSX repository on %osx_hostname% to SVN revision: %fds_revision%
plink %svn_logon% %scriptdir%/UPDATE_fds_onhost.csh  %fds_sourcedir% %fds_revision% %osx_hostname%
plink %svn_logon% %scriptdir%/UPDATE_fds_onhost_file.csh  %fds_builddir% makefile %fds_revision% %osx_hostname%

pause
