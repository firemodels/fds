@echo off

Rem Windows batch file to upload FDS test cases

set envfile=%homedrive%\%homepath%\fds_smv_env.bat
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

set windocs=docs_fds_smv_%docs_revision%.exe
set unixdocs=docs_fds_smv_%docs_revision%.tar.gz
echo Uploading %windocs%
echo Uploading %unixdocs%
pause

  set glabels=Type-Installer,Opsys-Windows,%docs_google_level%
  set dplatform=Windows
  set summary=Documentation  (SVN r%docs_revision% - Windows self-extracter)
  set file=%windocs%
  echo.
  echo Uploading %summary% 
       googlecode_upload.py --passwd-file-dir %google_password_dir% --config-dir none -s "%summary%" -p fds-smv -u gforney -l %glabels% %file%

echo.
echo Upload complete


  set glabels=Type-Archive,Opsys-Unix,%docs_google_level%
  set dplatform=Linux/OSX
  set summary=Documentation (SVN r%docs_revision% - gzipped/tar archive)
  set file=%unixdocs%
  echo.
  echo Uploading %summary% 
       googlecode_upload.py --passwd-file-dir %google_password_dir% --config-dir none -s "%summary%" -p fds-smv -u gforney -l %glabels% %file%

echo .
echo Upload complete

pause
