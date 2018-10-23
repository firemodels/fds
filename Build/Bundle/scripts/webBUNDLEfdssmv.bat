@echo off
setlocal EnableDelayedExpansion

set stopscript=0
call :getopts %*
if %stopscript% == 1 (
  exit /b
)
if NOT "%valid%" == "1" (
  call :usage
  exit /b
)

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
set env_defined=1
echo.
echo   Building FDS-Smokeview bundle for %platform%
Title  Building FDS-Smokeview bundle for %platform%

%svn_drive%

if "%platform%" == "windows" (
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

goto eof

:getopts
 if (%1)==() exit /b
 set valid=0
 set arg=%1
 if /I "%1" EQU "-help" (
   call :usage
   set stopscript=1
   exit /b
 )
 if /I "%1" EQU "windows" (
   set valid=1
 )
 if /I "%1" EQU "linux" (
   set valid=1
 )
 if /I "%1" EQU "osx" (
   set valid=1
 )
 shift
 if %valid% == 0 (
   echo.
   echo ***Error: the input argument %arg% is invalid
   echo.
   echo Usage:
   call :usage
   set stopscript=1
   exit /b
 )
if not (%1)==() goto getopts
exit /b

:usage  
echo webBUNDLEfdssmv [options] platform
echo. 
echo -help           - display this message
echo platform        - platform can be windows, linux or osx
exit /b

:eof
echo.
echo bundle build complete
pause
