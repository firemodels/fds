@echo off

Rem Windows batch file to upload Smokeview test files to
Rem the google download site.  This script assume that the Windows
Rem batch file, MAKEtest.bat, has already been run.

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
cd %svn_root%\SMV_5\for_bundle\to_google

Rem ----------------------------------------------------------
Rem should not need to edit any lines below

Rem set level=Release-1_Major
Rem set level=Release-2_Minor
set level=Release-3_Maintenance

set upload_win64=1
set upload_linux32=0
set upload_osx32=0

set version=%smv_version%

echo Uploading Smokeview %level% version=%smv_version% revision=%smv_revision%
pause

Rem --------------- 64 bit Windows ----------------
if not %upload_win64% == 1 goto endif_win64
  set glabels=Type-Installer,OpSys-Windows_64,%level%
  set dplatform=64 bit Windows
  set platform=win64
  set summary=Smokeview %smv_version% for %dplatform% (SVN r%smv_revision%)
  set exe=smv_%version%_%platform%.exe
  echo.
  echo Uploading %summary% - %exe%
  echo googlecode_upload.py --passwd-file-dir %google_password_dir% --config-dir none -s "%summary%" -p fds-smv -u %google_username% -l %glabels% %exe%
       googlecode_upload.py --passwd-file-dir %google_password_dir% --config-dir none -s "%summary%" -p fds-smv -u %google_username% -l %glabels% %exe%
:endif_win64

echo.
echo Uploads complete
pause
