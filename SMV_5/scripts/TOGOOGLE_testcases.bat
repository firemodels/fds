@echo off

Rem Windows batch file to upload FDS test cases

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
cd %svn_root%\Utilities\to_google

set winfile=verification_%verification_revision%.exe
set unixfile=verification_%verification_revision%.tar.gz

echo Uploading %winfile%
echo Uploading %unixfile%
pause

  set glabels=Type-Installer,Opsys-Windows,%verification_google_level%
  set dplatform=Windows
  set summary=Verification cases for FDS/Smokeview (SVN r%verification_revision% - Windows line endings)
  echo.
  echo Uploading %summary% - %wincase%
       googlecode_upload.py --passwd-file-dir %google_password_dir% --config-dir none -s "%summary%" -p fds-smv -u gforney -l %glabels% %winfile%

  set glabels=Type-Archive,Opsys-Unix,%verification_google_level%
  set dplatform=Linux/OSX
  set summary=Verification cases for FDS/Smokeview (SVN r%verification_revision% - Unix line endings)
  echo.
  echo Uploading %summary% - %wincase%
       googlecode_upload.py --passwd-file-dir %google_password_dir% --config-dir none -s "%summary%" -p fds-smv -u gforney -l %glabels% %unixfile%

echo.
echo Uploads complete
pause
