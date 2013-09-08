@echo off

Rem  Windows batch file to build Smokeview for all platforms.
Rem  This script builds LInux and OSX Smokeview's by doing a
Rem  remote shell (plink) to the NIST Linux cluster.

Rem setup environment variables (defining where repository resides etc) 

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

call %envfile%
echo Using the environment variables:
echo.
echo svn_root=%svn_root%
echo svn_drive=%svn_drive%
echo svn_logon=%svn_logon%
echo smv_version=%smv_version%
echo smv_revision=%smv_revision%
pause

Rem -----------------------------------------------------------
Rem shouldn't need to change any lines below

%svn_drive%
cd %svn_root%\smv\scripts
set version=%smv_version%_%smv_revision%

set scriptdir=FDS-SMV/SMV/scripts
set bundledir=FDS-SMV/SMV/for_bundle

plink %svn_logon% %scriptdir%/svn_update.csh
plink %svn_logon% %scriptdir%/make_smvs.csh %osx_hostname%
plink %svn_logon% %scriptdir%/make_dists.csh %version%

echo downloading Linux Smokeview files
pscp %svn_logon%:%bundledir%/smv_%version%_linux.tar.gz ..\for_bundle\uploads\.
pscp %svn_logon%:%bundledir%/smv_%version%_linux64.tar.gz ..\for_bundle\uploads\.

echo downloading MAC OSX Smokeview files
pscp %svn_logon%:%bundledir%/smv_%version%_osx.tar.gz ..\for_bundle\uploads\.

call make_smv_release_win32.bat %version%

echo uploading Windows Smokeview files
pscp  ..\for_bundle\smokeview_release.exe %svn_logon%:%bundledir%/.
pscp  ..\for_bundle\smokezip_release.exe %svn_logon%:%bundledir%/.
pause
