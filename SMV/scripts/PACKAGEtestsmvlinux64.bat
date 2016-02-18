@echo off
Title Building 64 bit linux test Smokeview bundle

Rem  Windows batch file to create an achive for a 64 bit Linux test smokeview

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

set scriptdir=%linux_svn_root%/SMV/scripts
set bundledir=%linux_svn_root%/SMV/uploads

echo making 64 bit Linux test distribution archive
plink %linux_logon% %scriptdir%/MAKEtestdistlinux64.csh %smv_revision% %linux_svn_root%

echo downloading Linux Smokeview files
pscp %linux_logon%:%bundledir%/smv_test_%smv_revision%_linux64.sh ..\uploads\.

pause
