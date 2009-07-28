@echo off

Rem Windows batch file to upload 64 bit OSX bundle to the google download site

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

call %envfile%

%svn_drive%
cd %svn_root%\Utilities\to_google

set glabels=Type-Installer,Opsys-OSX_64,%fds_google_level%
set dplatform=64 bit OSX
set summary=Bundled FDS and Smokeview for %dplatform% (SVN r%fds_revision%,%smv_revision%)
set exe=FDS_%fds_version%-SMV_%smv_version%_osx_64.tar.gz


echo Uploading %exe%
echo.
echo press any key to proceed with upload, CTRL c to abort
pause>NUL

  if not exist %exe% goto abort_upload
  echo.
  echo Uploading %summary% - %exe%
  echo googlecode_upload.py --passwd-file-dir %google_password_dir% --config-dir none -s "%summary%" -p fds-smv -u gforney -l %glabels% %exe%
       googlecode_upload.py --passwd-file-dir %google_password_dir% --config-dir none -s "%summary%" -p fds-smv -u gforney -l %glabels% %exe%

echo.
echo Uploads complete
pause
goto:eof

:abort_upload
echo %exe% does not exist - upload to google failed
pause
