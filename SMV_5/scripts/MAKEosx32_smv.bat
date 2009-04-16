@echo off

Rem  Windows batch file to build Smokeview for all platforms.
Rem  This script builds LInux and OSX Smokeview's by doing a
Rem  remote shell (plink) to the NIST Linux cluster.

Rem setup environment variables (defining where repository resides etc) 

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

call %envfile%
echo Using the environment variables:
echo.
echo Using SVN revision %smv_revision% to build Smokeview
pause

%svn_drive%
cd %svn_root%\smv_5\scripts
set version=%smv_version%_%smv_revision%

set scriptdir=FDS-SMV/SMV_5/scripts
set bundledir=FDS-SMV/SMV_5/for_bundle
set bindir=FDS-SMV/SMV_5/bin

plink %svn_logon% %scriptdir%/ssh_command.csh tiger.cfr.nist.gov %scriptdir% make_smv_osx.csh %smv_revision%

echo downloading compilation results for 32 bit MAC OSX Smokeview
pscp %svn_logon%:%bindir%/make_intel_osx_32.out %svn_root%\SMV_5\for_bundle\to_google\make_intel_osx_32.out

echo.
echo compilation complete
pause
