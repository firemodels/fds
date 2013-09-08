@echo off
Rem  Windows batch file to package a Linux Smokeview

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

%svn_drive%
cd %svn_root%\smv\scripts
set version=%smv_version%

set scriptdir=FDS-SMV/SMV/scripts
set bundledir=FDS-SMV/SMV/for_bundle/uploads

echo making Linux archives
plink %svn_logon% %scriptdir%/MAKEdistgen.csh %version% linux 32 %linux_hostname%

echo downloading Linux Smokeview files
pscp %svn_logon%:%bundledir%/smv_%version%_linux32.tar.gz ..\for_bundle\uploads\.

pause
