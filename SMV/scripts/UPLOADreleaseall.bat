@echo off

REM Windows batch file to upload Smokeview test files to
REM the download site.  This script assume that the Windows
REM batch file, MAKEtest.bat, has already been run.

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
cd %svn_root%\SMV\for_bundle\uploads

REM --------------- 32 bit Linux ----------------

  set platform=linux32
  set summary=Smokeview %smv_version% for %platform% (SVN r%smv_revision%)
  set exe=smv_%smv_version%_%platform%.tar.gz
  echo.
  echo Uploading %exe% 
       %upload% -k -ufds-smv:%api_key% -T %exe% https://api.bintray.com/content/%org_name%/%repo_name%/%package_name%/%smv_version%/%exe%;publish=1

echo.
echo Upload complete

REM --------------- 64 bit Linux ----------------

  set platform=linux64
  set summary=Smokeview %smv_version% for %platform% (SVN r%smv_revision%)
  set exe=smv_%smv_version%_%platform%.tar.gz
  echo.
  echo Uploading %exe% 
       %upload% -k -ufds-smv:%api_key% -T %exe% https://api.bintray.com/content/%org_name%/%repo_name%/%package_name%/%smv_version%/%exe%;publish=1

echo.
echo Upload complete

REM --------------- 32 bit OS X ----------------

  set platform=osx32
  set summary=Smokeview %smv_version% for %platform% (SVN r%smv_revision%)
  set exe=smv_%smv_version%_%platform%.tar.gz
  echo.
  echo Uploading %exe% 
       %upload% -k -ufds-smv:%api_key% -T %exe% https://api.bintray.com/content/%org_name%/%repo_name%/%package_name%/%smv_version%/%exe%;publish=1

echo.
echo Upload complete

REM --------------- 64 bit OS X ----------------

  set platform=osx64
  set summary=Smokeview %smv_version% for %platform% (SVN r%smv_revision%)
  set exe=smv_%smv_version%_%platform%.tar.gz
  echo.
  echo Uploading %exe% 
       %upload% -k -ufds-smv:%api_key% -T %exe% https://api.bintray.com/content/%org_name%/%repo_name%/%package_name%/%smv_version%/%exe%;publish=1

echo.
echo Upload complete

REM --------------- 32 bit Windows ----------------

set platform=win32
set summary=Smokeview %smv_version% for %platform% (SVN r%smv_revision%)
set exe=smv_%smv_version%_%platform%.exe
echo.
  echo Uploading %exe% 
     %upload% -k -ufds-smv:%api_key% -T %exe% https://api.bintray.com/content/%org_name%/%repo_name%/%package_name%/%smv_version%/%exe%;publish=1

echo.
echo Upload complete

REM --------------- 64 bit Windows ----------------

set platform=win64
set summary=Smokeview %smv_version% for %platform% (SVN r%smv_revision%)
set exe=smv_%smv_version%_%platform%.exe
echo.
  echo Uploading %exe% 
     %upload% -k -ufds-smv:%api_key% -T %exe% https://api.bintray.com/content/%org_name%/%repo_name%/%package_name%/%smv_version%/%exe%;publish=1

echo.
echo Upload complete
pause
