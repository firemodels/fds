@echo off

Rem Windows batch file to upload 32 bit Smokeview windows archive

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

  set level=Release-3_Maintenance
  set upload_win32=1
  set upload_linux32=0
  set upload_osx32=0

call %envfile%

echo Uploading Smokeview %level% version=%smv_version% revision=%smv_revision%
pause

%svn_drive%
cd %svn_root%\smv_5\for_bundle\to_google

Rem --------------- 32 bit Windows ----------------
if not %upload_win32% == 1 goto endif_win32
  set glabels=Type-Installer,OpSys-Windows,%level%
  set dplatform=32 bit Windows
  set platform=win32
  set summary=Smokeview %version% for %dplatform% (SVN r%smv_revision%)
  set exe=smv_%smv_version%_%platform%.exe
  echo.
  echo Uploading %summary% - %exe%
  echo googlecode_upload.py --passwd-file-dir %google_password_dir% --config-dir none -s "%summary%" -p fds-smv -u %google_username% -l %glabels% %exe%
       googlecode_upload.py --passwd-file-dir %google_password_dir% --config-dir none -s "%summary%" -p fds-smv -u %google_username% -l %glabels% %exe%
:endif_win32


echo.
echo Upload complete
pause
