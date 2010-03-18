@echo off

Rem Windows batch file to upload Smokeview osx files to
Rem the google download site.

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

set level=Release-3_Maintenance

set version=%smv_version%

echo Uploading Smokeview %level% version=%smv_version% revision=%smv_revision%
pause

Rem --------------- 32 bit OSX ----------------

  set glabels=Type-Archive,Opsys-OSX,%level%
  set dplatform=32 bit OSX
  set platform=osx32
  set summary=Smokeview %smv_version% for %dplatform% (SVN r%smv_revision%)
  set exe=smv_%smv_version%_%platform%.tar.gz
  echo.
  echo Uploading %summary% - %exe%
  echo googlecode_upload.py --passwd-file-dir %google_password_dir% --config-dir none  -s "%summary%" -p fds-smv -u %google_username% -l %glabels% %exe%
       googlecode_upload.py --passwd-file-dir %google_password_dir% --config-dir none  -s "%summary%" -p fds-smv -u %google_username% -l %glabels% %exe%

echo.
echo Uploads complete
pause
