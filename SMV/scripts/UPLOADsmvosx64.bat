@echo off

:: Windows batch file to upload Smokeview test files to
:: the download site.  This script assume that the Windows
:: batch file, MAKEtest.bat, has already been run.

set platform=osx64

set envfile=%userprofile%\fds_smv_env.bat
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

set summary=Smokeview %smv_version% for %platform% (GIT r%smv_revision%)
set exe=smv_%smv_version%_%platform%.sh

echo Uploading %exe%
echo.
echo Press any key to proceed with upload, CTRL c to abort
pause>NUL

  if not exist %exe% goto abort_upload
  echo.
  echo Uploading %summary% - %exe%

  if not exist %bintray_api_key% goto no_key_file
  echo.
  set /p api_key=<%bintray_api_key%

:: Repository information

  set org_name=nist-fire-research
  set repo_name=releases
  set package_name=SMV

:: Upload and publish file

  %upload% -k -ufds-smv:%api_key% -T %exe% https://api.bintray.com/content/%org_name%/%repo_name%/%package_name%/%smv_version%/%exe%;publish=1

echo.
echo Upload complete
pause
goto:eof

:abort_upload
echo Error: %exe% does not exist - upload failed.
pause
goto:eof

:no_key_file
echo Error: Bintray API key does not exist in %bintray_api_key% - upload failed.
pause
goto:eof
