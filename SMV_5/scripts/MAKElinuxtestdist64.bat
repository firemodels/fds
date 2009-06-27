@echo off

Rem  Windows batch file to create an achive for a 64 bit Linux test smokeview

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

%svn_drive%
cd %svn_root%\smv_5\scripts
set version=%smv_version%

set scriptdir=FDS-SMV/SMV_5/scripts
set bundledir=FDS-SMV/SMV_5/for_bundle

echo making 64 bit Linux test distribution archive
plink %svn_logon% %scriptdir%/make_linux64_test_dist.csh %version% %smv_revision%

echo downloading Linux Smokeview files
pscp %svn_logon%:%bundledir%/smv_%version%_%smv_revision%_linux64.tar.gz ..\for_bundle\to_google\.
pause
