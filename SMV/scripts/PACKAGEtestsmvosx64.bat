@echo off
Title Packaging test Smokeview for 64 bit OSX

Rem  Windows batch file to create an OSX achive for an OSX test smokeview

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

cd "%svn_root%\..\Google Drive\SMV_Test_Versions
set gupload=%CD%

cd %svn_root%\smv\scripts

set scriptdir=%linux_svn_root%/SMV/scripts
set bundledir=%linux_svn_root%/SMV/uploads

echo making OSX test distribution archive
plink %svn_logon% %scriptdir%/MAKEtestdistosx64.csh %smv_revision% %osx_hostname% %linux_svn_root%

echo downloading OSX test distribution archive
pscp %svn_logon%:%bundledir%/smv_test_%smv_revision%_osx64.sh ..\uploads\.

echo copying ..\uploads\smv_test_%smv_revision%_osx64.sh to %gupload%
copy ..\uploads\smv_test_%smv_revision%_osx64.sh "%gupload%"

pause
