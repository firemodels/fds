@echo off

Rem  Windows batch file to build Smokeview for all platforms.
Rem  This script builds LInux and OSX Smokeview's by doing a
Rem  remote shell (plink) to the NIST Linux cluster.

Rem setup environment variables (defining where repository resides etc) 

set envfile=%homedrive%\%homepath%\fds_smv_env.bat
IF EXIST %envfile% GOTO endif_envexist
echo ***Fatal error.  The environment setup file %envfile% does not exist. 
echo Create a file named %envfile% and use SMV_5/scripts/fds_smv_env_template.bat
echo as an example.
echo.
echo Aborting now...
pause>NUL
goto:eof

:endif_envexist

call %envfile%

%svn_drive%
cd %svn_root%\smv_5\scripts
set version=%smv_version%_%smv_revision%

set scriptdir=FDS-SMV/SMV_5/scripts
set bundledir=FDS-SMV/SMV_5/for_bundle

echo making Linux and OSX distribution archives
plink %svn_logon% %scriptdir%/make_dists.csh %version%

echo downloading Linux Smokeview files
pscp %svn_logon%:%bundledir%/smv_%version%_linux.tar.gz ..\for_bundle\to_google\.
Rem pscp %svn_logon%:%bundledir%/smv_%version%_linux_64.tar.gz ..\for_bundle\to_google\.
Rem pscp %svn_logon%:%scriptdir%/make_intel_linux_32.out ..\for_bundle\to_google\.
Rem pscp %svn_logon%:%scriptdir%/make_intel_linux_64.out ..\for_bundle\to_google\.

echo downloading MAC OSX Smokeview files
pscp %svn_logon%:%bundledir%/smv_%version%_osx.tar.gz ..\for_bundle\to_google\.
Rem pscp %svn_logon%:%scriptdir%/make_intel_osx_32.out ..\for_bundle\to_google\.

pause
