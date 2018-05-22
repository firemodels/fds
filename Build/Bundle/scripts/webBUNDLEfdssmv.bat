@echo off
setlocal EnableDelayedExpansion
set platform=%1

:: batch file to generate Windows, Linux or OSX FDS-SMV bundles

:: setup environment variables (defining where repository resides etc) 

set envfile="%userprofile%"\fds_smv_env.bat
IF EXIST %envfile% GOTO endif_envexist
echo ***Fatal error.  The environment setup file %envfile% does not exist. 
echo Create a file named %envfile% and use smv/scripts/fds_smv_env_template.bat
echo as an example.
echo.
echo Aborting now...
pause>NUL
goto:eof

:endif_envexist

call %envfile%
echo.
echo   Building FDS-Smokeview bundle for %platform%
Title  Building FDS-Smokeview bundle for %platform%

%svn_drive%

if "%platform%" == "windows" (
  set platform=64
  call "%svn_root%\fds\Build\Bundle\windows\make_bundle"
  goto eof
)
if "%platform%" == "linux" (
  set bundledir=%fds_version%-%smv_version%_linux64
  plink %linux_logon% %linux_svn_root%/fds/Build/Bundle/linux/make_bundle_fromweb.sh

  echo Downloading compressed archive to:
  echo   %svn_root%\fds\Build\Bundle\uploads\!bundledir!.sh
  pscp %linux_logon%:%linux_svn_root%/fds/Build/Bundle/uploads/!bundledir!.sh   %svn_root%\fds\Build\Bundle\uploads\.
  pscp %linux_logon%:%linux_svn_root%/fds/Build/Bundle/uploads/!bundledir!.sha1 %svn_root%\fds\Build\Bundle\uploads\.
  goto eof
)
if "%platform%" == "osx" (
  set bundledir=%fds_version%-%smv_version%_osx64
  plink %osx_logon% %linux_svn_root%/fds/Build/Bundle/osx/make_bundle_fromweb.sh

  echo Downloading compressed archive to:
  echo   %svn_root%\fds\Build\Bundle\uploads\!bundledir!.sh
  pscp %osx_logon%:%linux_svn_root%/fds/Build/Bundle/uploads/!bundledir!.sh   %svn_root%\fds\Build\Bundle\uploads\.
  pscp %osx_logon%:%linux_svn_root%/fds/Build/Bundle/uploads/!bundledir!.sha1 %svn_root%\fds\Build\Bundle\uploads\.
  goto eof
)

:eof
echo.
echo bundle build complete
pause
