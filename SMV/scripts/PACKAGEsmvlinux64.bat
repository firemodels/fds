@echo off
Title Building 64 bit Linux release Smokeview bundle

Rem  Windows batch file to create an achive for a 64 bit Linux smokeview

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

set scriptdir=%linux_svn_root%/SMV/scripts
set bundledir=%linux_svn_root%/SMV/uploads

echo.
echo --- making 64 bit Linux Smokeview installer ---
echo.
plink %linux_logon% %scriptdir%/MAKEdistgen.sh %version% linux 64 %linux_hostname% %fds_edition% %linux_svn_root%

echo.
echo --- downloading Linux Smokeview installer ---
echo.
pscp %linux_logon%:%bundledir%/smv_%version%_linux64.sh ..\uploads\.

pause
