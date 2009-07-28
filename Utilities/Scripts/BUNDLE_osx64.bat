@echo off
Rem Batch file used to bundle a 64 bit OSX of FDS/Smokeview

set envfile="%homedrive%\%homepath%\fds_smv_env.bat"
IF EXIST %envfile% GOTO endif_envexist
echo ***Fatal error.  The environment setup file %envfile% does not exist. 
echo Create a file named %envfile% and use SMV_5/scripts/fds_smv_env_template.bat
echo as an example.
echo.
echo Aborting now...
pause>NUL
goto:eof

:endif_envexist

Rem location of batch files used to set up Intel compilation environment

call %envfile%

set bundledir=FDS_%fds_version%-SMV_%smv_version%_osx_64
plink %svn_logon% %linux_svn_root%/Utilities/Scripts/bundle_osx_64.csh %linux_svn_root% %bundledir%

echo Downloading compressed archive to:
echo   %svn_root%\Utilities\to_google\%bundledir%.tar.gz
pscp %svn_logon%:%linux_svn_root%/Utilities/to_google/%bundledir%.tar.gz %svn_root%/Utilities/to_google/.

pause