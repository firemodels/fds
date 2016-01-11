@echo off

Rem  Windows batch file to create an OSX achive for an OSX smokeview

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
set version=%smv_version%

set scriptdir=%linux_svn_root%/SMV/scripts
set bundledir=%linux_svn_root%/SMV/uploads

echo making 64 bit Smokeview OSX distribution archive
plink %svn_logon% %scriptdir%/MAKEdistgen.csh %version% osx 64 %osx_hostname% %fds_edition% %linux_svn_root%

echo downloading 64 bit Smokeview OSX distribution archive
pscp %svn_logon%:%bundledir%/smv_%version%_osx64.sh ..\uploads\.

echo copying ..\uploads\smv_%version%_osx64.sh to %gupload%
copy ..\uploads\smv_%version%_osx64.sh "%gupload%"

pause
