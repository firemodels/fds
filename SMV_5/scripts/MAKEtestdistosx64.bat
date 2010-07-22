@echo off

Rem  Windows batch file to create an OSX achive for an OSX test smokeview

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

set scriptdir=FDS-SMV/SMV_5/scripts
set bundledir=FDS-SMV/SMV_5/for_bundle

echo making OSX test distribution archive
plink %svn_logon% %scriptdir%/MAKEtestdistosx64.csh %smv_revision% %OSXHOST%

echo downloading OSX test distribution archive
pscp %svn_logon%:%bundledir%/smv_test_%smv_revision%_osx_64.tar.gz ..\for_bundle\to_google\.


pause
