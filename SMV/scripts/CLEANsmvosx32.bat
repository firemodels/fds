@echo off
Title Cleaning Smokeview for 32 bit OSX

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

set scriptdir=%linux_svn_root%/SMV/scripts
set smvdir=%linux_svn_root%/SMV/Build/

plink %svn_logon% %scriptdir%/CLEAN_smv_onhost.csh %smvdir%/INTEL_OSX_32 %osx_hostname% clean

pause